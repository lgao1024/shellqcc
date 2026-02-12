@if (@X)==(@Y) @end /*
@echo off
setlocal EnableExtensions

if "%~1"=="" (
  echo Usage: %~nx0 input_file [output_file]
  exit /b 2
)

where cscript >nul 2>nul
if errorlevel 1 (
  echo ERROR: cscript not found. Windows Script Host may be disabled.
  exit /b 2
)

cscript //nologo //E:JScript "%~f0" "%~f1" "%~2"
exit /b %ERRORLEVEL%
*/

//
// =============== JScript section (runs under cscript) ===============
// Usage: cscript //nologo //E:JScript xtal_convert.bat input [output]
//
function die(msg){ WScript.StdErr.WriteLine(msg); WScript.Quit(2); }

if (WScript.Arguments.length < 1) die("Usage: xtal_convert.bat input_file [output_file]");

var inPath  = WScript.Arguments.Item(0);
var outPath = (WScript.Arguments.length >= 2) ? WScript.Arguments.Item(1) : "";

var fso = new ActiveXObject("Scripting.FileSystemObject");
if (!fso.FileExists(inPath)) die("ERROR: input not found: " + inPath);

function extLower(p){
  var i = p.lastIndexOf(".");
  if (i < 0) return "";
  return p.substring(i).toLowerCase();
}
function changeExt(p, newExt){
  var i = p.lastIndexOf(".");
  if (i < 0) return p + newExt;
  return p.substring(0, i) + newExt;
}

// -------- classification for sorting --------
// Nonmetals + noble gases + metalloids treated as nonmetal group
var NONMETAL = {
  "H":1,"C":1,"N":1,"O":1,"F":1,"P":1,"S":1,"Se":1,
  "Cl":1,"Br":1,"I":1,"At":1,"Ts":1,
  "He":1,"Ne":1,"Ar":1,"Kr":1,"Xe":1,"Rn":1,"Og":1,
  "B":1,"Si":1,"Ge":1,"As":1,"Sb":1,"Te":1,"Po":1
};
function classKey(elem){ return NONMETAL.hasOwnProperty(elem) ? 1 : 0; } // metal=0 first

// -------- helpers --------
function trim(s){ return String(s).replace(/^\s+|\s+$/g, ""); }
function stripComment(line){ return trim(String(line).replace(/#.*$/, "")); }
function splitWS(line){
  var t = trim(line);
  if (!t) return [];
  return t.split(/\s+/);
}
function unquote(s){
  s = trim(s);
  if (s.length >= 2) {
    var a = s.charAt(0), b = s.charAt(s.length-1);
    if ((a=="'" && b=="'") || (a=='"' && b=='"')) return s.substring(1, s.length-1);
  }
  return s;
}
function parseNumber(s){
  if (s === null || s === undefined) throw "empty";
  var t = trim(String(s));
  if (!t) throw "empty";

  // 0.123(4) -> 0.123
  t = t.replace(/\(\d+\)$/, "");

  // fraction a/b
  if (/^[+-]?\d+\/\d+$/.test(t)) {
    var p = t.split("/");
    var n = parseFloat(p[0]);
    var d = parseFloat(p[1]);
    if (d === 0) throw "div0";
    return n / d;
  }

  var v = parseFloat(t);
  if (isNaN(v)) throw "nan";
  return v;
}
function isNumberLike(s){
  try { parseNumber(s); return true; } catch(e){ return false; }
}
function fmt6(x){
  var v = Number(x);
  return v.toFixed(6);
}
function parseTxtHeader(line){
  line = trim(line);
  var sg="", rest="";
  if (line.charAt(0)=="'") {
    var j = line.indexOf("'", 1);
    if (j < 0) throw "bad header";
    sg = line.substring(1, j);
    rest = trim(line.substring(j+1));
  } else if (line.charAt(0)=='"') {
    var k = line.indexOf('"', 1);
    if (k < 0) throw "bad header";
    sg = line.substring(1, k);
    rest = trim(line.substring(k+1));
  } else {
    var p = splitWS(line);
    if (p.length < 7) throw "bad header";
    sg = p[0];
    rest = p.slice(1).join(" ");
  }

  var a = splitWS(rest);
  if (a.length < 6) throw "bad header";
  return {
    sg: sg,
    a: parseNumber(a[0]), b: parseNumber(a[1]), c: parseNumber(a[2]),
    alpha: parseNumber(a[3]), beta: parseNumber(a[4]), gamma: parseNumber(a[5])
  };
}

function sortAtoms(atoms){
  atoms.sort(function(u, v){
    if (u.classKey !== v.classKey) return u.classKey - v.classKey; // metal first
    if (u.elem < v.elem) return -1;
    if (u.elem > v.elem) return 1;
    if (u.label < v.label) return -1;
    if (u.label > v.label) return 1;
    return 0;
  });
}

function readAllLines(path){
  var ts = fso.OpenTextFile(path, 1, false); // ForReading
  var out = [];
  while (!ts.AtEndOfStream) out.push(ts.ReadLine());
  ts.Close();
  return out;
}
function writeText(path, content){
  var ts = fso.OpenTextFile(path, 2, true); // ForWriting, create
  ts.Write(content);
  ts.Close();
}

// -------- TXT -> CIF --------
function txt2cif(inFile, outFile){
  var lines = readAllLines(inFile);
  var header = null;
  var atoms = [];
  var elemCount = {};

  for (var i=0; i<lines.length; i++){
    var line = stripComment(lines[i]);
    if (!line) continue;

    if (!header){
      try { header = parseTxtHeader(line); }
      catch(e){ die("ERROR: first line must be: sg a b c alpha beta gamma. Got: " + line); }
      continue;
    }

    var tok = splitWS(line);
    if (tok.length < 4) die("ERROR: atom line has too few columns: " + line);

    var label="", elem="", x=0, y=0, z=0, occ=1.0;

    if (isNumberLike(tok[1])) {
      // element x y z [occ]
      elem = tok[0];
      x = parseNumber(tok[1]); y = parseNumber(tok[2]); z = parseNumber(tok[3]);
      occ = (tok.length >= 5) ? parseNumber(tok[4]) : 1.0;
      if (!elemCount.hasOwnProperty(elem)) elemCount[elem]=0;
      elemCount[elem] += 1;
      label = elem + elemCount[elem];
    } else {
      // label element x y z [occ]
      if (tok.length < 5) die("ERROR: expected at least 5 fields: label element x y z [occ]. Got: " + line);
      label = tok[0]; elem = tok[1];
      x = parseNumber(tok[2]); y = parseNumber(tok[3]); z = parseNumber(tok[4]);
      occ = (tok.length >= 6) ? parseNumber(tok[5]) : 1.0;
    }

    atoms.push({ label:label, elem:elem, x:x, y:y, z:z, occ:occ, classKey:classKey(elem) });
  }

  if (!header) die("ERROR: no valid first line found in " + inFile);

  sortAtoms(atoms);

  var dataName = fso.GetBaseName(outFile).replace(/\s+/g, "_");
  var out = "";
  out += "data_" + dataName + "\r\n\r\n";
  out += "_symmetry_space_group_name_H-M  '" + header.sg + "'\r\n\r\n";
  out += "_cell_length_a    " + fmt6(header.a) + "\r\n";
  out += "_cell_length_b    " + fmt6(header.b) + "\r\n";
  out += "_cell_length_c    " + fmt6(header.c) + "\r\n";
  out += "_cell_angle_alpha " + fmt6(header.alpha) + "\r\n";
  out += "_cell_angle_beta  " + fmt6(header.beta) + "\r\n";
  out += "_cell_angle_gamma " + fmt6(header.gamma) + "\r\n\r\n";
  out += "loop_\r\n";
  out += "_atom_site_label\r\n";
  out += "_atom_site_type_symbol\r\n";
  out += "_atom_site_fract_x\r\n";
  out += "_atom_site_fract_y\r\n";
  out += "_atom_site_fract_z\r\n";
  out += "_atom_site_occupancy\r\n";

  for (var j=0; j<atoms.length; j++){
    var a = atoms[j];
    out += a.label + " " + a.elem + " " + fmt6(a.x) + " " + fmt6(a.y) + " " + fmt6(a.z) + " " + fmt6(a.occ) + "\r\n";
  }
  out += "\r\n";

  writeText(outFile, out);
  WScript.Echo("Wrote: " + outFile);
}

// -------- CIF -> TXT --------
function cif2txt(inFile, outFile){
  var lines = readAllLines(inFile);
  var sg = "";
  var cell = {a:null,b:null,c:null,alpha:null,beta:null,gamma:null};

  var inLoop=false, reading=false;
  var cols=[], idx={label:-1,type:-1,x:-1,y:-1,z:-1,occ:-1};
  var atoms=[];
  var elemCount={};

  function setCell(tag, val){
    try { cell[tag] = parseNumber(val); } catch(e){}
  }

  for (var i=0; i<lines.length; i++){
    var line = stripComment(lines[i]);
    if (!line) continue;
    var low = line.toLowerCase();

    if (low.indexOf("_symmetry_space_group_name_h-m")===0) {
      sg = unquote(line.replace(/^[^\s]+\s+/, ""));
      continue;
    }
    if (!sg && low.indexOf("_space_group_name_h-m_alt")===0) {
      sg = unquote(line.replace(/^[^\s]+\s+/, ""));
      continue;
    }

    if (low.indexOf("_cell_length_a")===0) { setCell("a", line.replace(/^[^\s]+\s+/, "")); continue; }
    if (low.indexOf("_cell_length_b")===0) { setCell("b", line.replace(/^[^\s]+\s+/, "")); continue; }
    if (low.indexOf("_cell_length_c")===0) { setCell("c", line.replace(/^[^\s]+\s+/, "")); continue; }
    if (low.indexOf("_cell_angle_alpha")===0) { setCell("alpha", line.replace(/^[^\s]+\s+/, "")); continue; }
    if (low.indexOf("_cell_angle_beta")===0) { setCell("beta", line.replace(/^[^\s]+\s+/, "")); continue; }
    if (low.indexOf("_cell_angle_gamma")===0) { setCell("gamma", line.replace(/^[^\s]+\s+/, "")); continue; }

    if (low === "loop_") {
      inLoop=true; reading=false;
      cols=[]; idx={label:-1,type:-1,x:-1,y:-1,z:-1,occ:-1};
      continue;
    }

    if (inLoop && !reading && line.charAt(0)=="_") {
      var c = low;
      cols.push(c);
      var pos = cols.length-1;
      if (c==="_atom_site_label") idx.label=pos;
      if (c==="_atom_site_type_symbol") idx.type=pos;
      if (c==="_atom_site_fract_x") idx.x=pos;
      if (c==="_atom_site_fract_y") idx.y=pos;
      if (c==="_atom_site_fract_z") idx.z=pos;
      if (c==="_atom_site_occupancy") idx.occ=pos;
      continue;
    }

    if (inLoop && !reading && line.charAt(0)!="_") {
      if (idx.x>=0 && idx.y>=0 && idx.z>=0) reading=true;
      else { inLoop=false; reading=false; continue; }
      // fall through to parse this as data
    }

    if (inLoop && reading) {
      if (low==="loop_" || line.charAt(0)=="_" || low.indexOf("data_")===0) { inLoop=false; reading=false; continue; }

      var f = splitWS(line);
      if (f.length < cols.length) continue;

      var label = (idx.label>=0) ? f[idx.label] : "";
      var elem  = (idx.type>=0)  ? f[idx.type]  : "";
      if (!elem && label) elem = label.replace(/\d+$/, "");
      if (!elem) elem = "X";

      var x = parseNumber(f[idx.x]), y = parseNumber(f[idx.y]), z = parseNumber(f[idx.z]);
      var occ = 1.0;
      if (idx.occ>=0) {
        var o = f[idx.occ];
        if (o !== "." && o !== "?") occ = parseNumber(o);
      }

      if (!label) {
        if (!elemCount.hasOwnProperty(elem)) elemCount[elem]=0;
        elemCount[elem] += 1;
        label = elem + elemCount[elem];
      }

      atoms.push({label:label, elem:elem, x:x, y:y, z:z, occ:occ, classKey:classKey(elem)});
    }
  }

  if (!sg) sg="P1";
  if (cell.a===null) cell.a=1; if (cell.b===null) cell.b=1; if (cell.c===null) cell.c=1;
  if (cell.alpha===null) cell.alpha=90; if (cell.beta===null) cell.beta=90; if (cell.gamma===null) cell.gamma=90;

  sortAtoms(atoms);

  var sgOut = (/\s/.test(sg) ? ("'" + sg + "'") : sg);
  var out = "";
  out += sgOut + " " + fmt6(cell.a) + " " + fmt6(cell.b) + " " + fmt6(cell.c) + " " + fmt6(cell.alpha) + " " + fmt6(cell.beta) + " " + fmt6(cell.gamma) + "\r\n";
  for (var j=0; j<atoms.length; j++){
    var a = atoms[j];
    out += a.label + " " + a.elem + " " + fmt6(a.x) + " " + fmt6(a.y) + " " + fmt6(a.z) + " " + fmt6(a.occ) + "\r\n";
  }
  writeText(outFile, out);
  WScript.Echo("Wrote: " + outFile);
}

// -------- dispatch by suffix --------
var ext = extLower(inPath);
if (!outPath) {
  outPath = (ext === ".cif") ? changeExt(inPath, ".txt") : changeExt(inPath, ".cif");
}
if (ext === ".cif") cif2txt(inPath, outPath);
else txt2cif(inPath, outPath);
