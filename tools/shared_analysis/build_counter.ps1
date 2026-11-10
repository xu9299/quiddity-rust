param(
    [Parameter(Mandatory=$true)][string]$DynamoRIO,
    [Parameter(Mandatory=$true)][string]$OutputDirectory
)
$ErrorActionPreference = 'Stop'
$DynamoRIO = (Resolve-Path -LiteralPath $DynamoRIO).Path
$OutputDirectory = [IO.Path]::GetFullPath($OutputDirectory)
New-Item -ItemType Directory -Force -Path $OutputDirectory | Out-Null
# Run in an x64 MSVC developer environment, with its Windows SDK include paths.
& cl.exe /nologo /c /W4 /WX /O2 /GS- /TC /DWINDOWS /DX86_64 `
    "/I$DynamoRIO/include" "/I$DynamoRIO/ext/include" `
    "/Fo$OutputDirectory/instruction_counter.obj" "$PSScriptRoot/instruction_counter.c"
if ($LASTEXITCODE -ne 0) { throw 'Counter compilation failed' }
& link.exe /nologo /DLL /NOENTRY /NODEFAULTLIB /MANIFEST:NO /OPT:REF /OPT:ICF `
    "/OUT:$OutputDirectory/instruction_counter.dll" "$OutputDirectory/instruction_counter.obj" `
    "$DynamoRIO/lib64/release/dynamorio.lib" "$DynamoRIO/ext/lib64/release/drmgr.lib" `
    "$DynamoRIO/ext/lib64/release/drwrap.lib"
if ($LASTEXITCODE -ne 0) { throw 'Counter link failed' }
& ml64.exe /nologo /c "/Fo$OutputDirectory/calibration.obj" "$PSScriptRoot/calibration.asm"
if ($LASTEXITCODE -ne 0) { throw 'Calibration assembly failed' }
& link.exe /nologo /DLL /NOENTRY /NODEFAULTLIB /MANIFEST:NO /EXPORT:icount_loop /EXPORT:icount_marker `
    "/OUT:$OutputDirectory/icount_calibration.dll" "$OutputDirectory/calibration.obj"
if ($LASTEXITCODE -ne 0) { throw 'Calibration link failed' }
