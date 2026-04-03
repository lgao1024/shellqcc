@if (@CodeSection == @Batch) @then
@echo off
setlocal

:: This is a Hybrid Batch + JScript file.
:: It uses cscript to run the JScript engine on this exact same file.
cscript //e:JScript //nologo "%~f0" %*

goto :EOF
@end

// ==========================================
// JScript Code Starts Here
// ==========================================
var args = WScript.Arguments;
var fso = new ActiveXObject("Scripting.FileSystemObject");

// ==========================================
// Help Information and Input Validation
// ==========================================
if (args.length === 0) {
    WScript.Echo("Usage:");
    WScript.Echo("  1. Split file: cif_manager.bat <merged_cif_file.cif>");
    WScript.Echo("  2. Merge files: cif_manager.bat <file1.cif> <file2.cif> ... <fileN.cif>");
    WScript.Quit(1);
}

// ==========================================
// Function 1: Split a single merged CIF file
// ==========================================
if (args.length === 1) {
    var inputFile = args(0);
    
    if (!fso.FileExists(inputFile)) {
        WScript.Echo("Error: File '" + inputFile + "' not found");
        WScript.Quit(1);
    }

    var baseName = fso.GetBaseName(inputFile);
    var outDir = baseName + "_splited";
    
    if (!fso.FolderExists(outDir)) {
        fso.CreateFolder(outDir);
    }

    WScript.Echo("Starting to split CIF files into directory: " + outDir + "\\");

    var file = fso.OpenTextFile(inputFile, 1, false); // 1 = ForReading
    var counter = 0;
    var hasData = false;
    var mainBuffer = [];
    var tempBuffer = [];
    var ccdcNum = "";

    var rxData = /^data_/;
    var rxCommentOrEmpty = /^(\s*#|\s*$)/;
    var rxCcdc = /_database_code_depnum_ccdc_archive/;
    var rxDigits = /[0-9]+/;

    function flushBlock() {
        if (hasData) {
            var outName = (ccdcNum !== "") ? ccdcNum : ("structure" + counter);
            var outPath = fso.BuildPath(outDir, outName + ".cif");
            var outFile = fso.CreateTextFile(outPath, true); // true = overwrite
            outFile.Write(mainBuffer.join("\n"));
            outFile.Close();
            WScript.Echo("Generated: " + outPath);
        }
    }

    while (!file.AtEndOfStream) {
        var line = file.ReadLine();

        // When a new data_ tag is matched
        if (rxData.test(line)) {
            flushBlock();
            counter++;
            ccdcNum = "";
            hasData = true;
            
            mainBuffer = [];
            if (tempBuffer.length > 0) {
                mainBuffer = mainBuffer.concat(tempBuffer);
            }
            mainBuffer.push(line);
            tempBuffer = [];
            continue;
        }

        // Process global comments at the very beginning of the file
        if (!hasData) {
            tempBuffer.push(line);
            continue;
        }

        // If it is a comment line or a blank line, put it into temp buffer
        if (rxCommentOrEmpty.test(line)) {
            tempBuffer.push(line);
        } else {
            // Normal CIF data line, flush temp buffer to main buffer
            if (tempBuffer.length > 0) {
                mainBuffer = mainBuffer.concat(tempBuffer);
                tempBuffer = [];
            }
            mainBuffer.push(line);
            
            // Extract the CCDC number
            if (rxCcdc.test(line)) {
                var match = rxDigits.exec(line);
                if (match) {
                    ccdcNum = match[0];
                }
            }
        }
    }

    // Process the last data block
    if (hasData) {
        if (tempBuffer.length > 0) {
            mainBuffer = mainBuffer.concat(tempBuffer);
        }
        flushBlock();
    }

    file.Close();
    WScript.Echo("Split complete!");
    WScript.Quit(0);
}

// ==========================================
// Function 2: Merge multiple CIF files
// ==========================================
var files = [];
for (var i = 0; i < args.length; i++) {
    if (!fso.FileExists(args(i))) {
        WScript.Echo("Error: File '" + args(i) + "' not found");
        WScript.Quit(1);
    }
    files.push(args(i));
}

var currentOrder = [];
for (var i = 0; i < files.length; i++) {
    currentOrder.push(i + 1);
}

function parseOrder(inputStr) {
    var newOrder = [];
    var parts = inputStr.replace(/,/g, ' ').split(/\s+/);
    for (var i = 0; i < parts.length; i++) {
        var p = parts[i];
        if (!p) continue;
        
        // Match hyphen range, e.g., 5-6 or 3-1
        if (/^[0-9]+-[0-9]+$/.test(p)) {
            var bounds = p.split('-');
            var start = parseInt(bounds[0], 10);
            var end = parseInt(bounds[1], 10);
            if (start <= end) {
                for (var j = start; j <= end; j++) newOrder.push(j);
            } else {
                for (var j = start; j >= end; j--) newOrder.push(j);
            }
        } 
        // Match pure numbers
        else if (/^[0-9]+$/.test(p)) {
            newOrder.push(parseInt(p, 10));
        }
    }
    return newOrder;
}

while (true) {
    WScript.Echo("\nCurrent files to be merged and their order:");
    for (var i = 0; i < currentOrder.length; i++) {
        var idx = currentOrder[i] - 1;
        WScript.Echo("  " + (i + 1) + ": " + files[idx]);
    }

    WScript.Echo("\nYou can enter the desired merge order (e.g., 1, 3, 5-6, 2).");
    WScript.StdOut.Write("Press Enter to keep the current order and merge directly: ");
    
    var userInput = WScript.StdIn.ReadLine();

    if (!userInput || userInput.replace(/^\s+|\s+$/g, '') === "") {
        break; // Empty input, keep current order
    }

    var tempOrder = parseOrder(userInput);
    var valid = true;
    for (var i = 0; i < tempOrder.length; i++) {
        var idx = tempOrder[i];
        if (idx < 1 || idx > files.length) {
            WScript.Echo("Warning: Index " + idx + " is out of range! Please try again.");
            valid = false;
            break;
        }
    }

    if (valid && tempOrder.length > 0) {
        WScript.Echo("\n---- Preview of updated file order ----");
        for (var i = 0; i < tempOrder.length; i++) {
            var idx = tempOrder[i] - 1;
            WScript.Echo("  " + (i + 1) + ": " + files[idx]);
        }
        
        WScript.StdOut.Write("Confirm to merge with this order? (y/n, Enter to confirm): ");
        var confirm = WScript.StdIn.ReadLine();
        if (!confirm || /^[Yy]/.test(confirm)) {
            currentOrder = tempOrder;
            break;
        }
    } else if (tempOrder.length === 0) {
        WScript.Echo("Warning: Invalid input format.");
    }
}

// Generate output filename with timestamp
function pad(n) { return n < 10 ? '0' + n : n; }
var d = new Date();
var timestamp = d.getFullYear() + pad(d.getMonth() + 1) + pad(d.getDate()) + "_" + 
                pad(d.getHours()) + pad(d.getMinutes()) + pad(d.getSeconds());
var mergedOutput = "merged_output_" + timestamp + ".cif";

WScript.Echo("\nStarting merge...");
var outFile = fso.CreateTextFile(mergedOutput, true);

for (var i = 0; i < currentOrder.length; i++) {
    var idx = currentOrder[i] - 1;
    var inFile = fso.OpenTextFile(files[idx], 1, false);
    
    if (!inFile.AtEndOfStream) {
        var content = inFile.ReadAll();
        outFile.Write(content);
    }
    inFile.Close();
    
    // Ensure newline separation
    outFile.Write("\n\n");
}

outFile.Close();
WScript.Echo("Merge complete! Successfully output to file: " + mergedOutput);