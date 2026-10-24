# Render the concept figure's 3D panels with headless Edge.
# Usage (repository root, explorer data built and served on port 8000 by scripts/explorer.sh):
#   powershell -File scripts/figure_renders/render.ps1 [-Id 5SAK_ZRY] [-Out docs/figures/renders]
# Writes one transparent PNG per scene: input, receptor, box, engine, evaluation.
param(
    [string]$Id = "5SAK_ZRY",
    [string]$Out = "docs/figures/renders",
    [string]$Edge = "C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
)
# Not "Stop": Edge logs harmless messages to stderr, which PowerShell 5.1 would turn into errors.
$ErrorActionPreference = "Continue"
$page = "explorer/zz_render.html"
Copy-Item "scripts/figure_renders/render.html" $page
New-Item -ItemType Directory -Force $Out | Out-Null
$profile = Join-Path $env:TEMP "dockprep-render-profile"
try {
    foreach ($scene in "input", "receptor", "box", "engine", "evaluation") {
        # A fresh profile per scene: reusing one left later scenes' canvases blank.
        $dom = & $Edge --headless=new --user-data-dir="$profile-$scene" --virtual-time-budget=20000 `
            --window-size=900,560 --dump-dom "http://localhost:8000/zz_render.html?id=$Id&scene=$scene" 2>$null |
            Out-String
        if ($dom -notmatch 'data:image/png;base64,([A-Za-z0-9+/=\s]+)<') { throw "No image for scene $scene" }
        $bytes = [Convert]::FromBase64String(($Matches[1] -replace '\s', ''))
        [IO.File]::WriteAllBytes((Join-Path (Resolve-Path $Out) "$scene.png"), $bytes)
        Write-Output "$scene.png $($bytes.Length) bytes"
    }
} finally {
    Remove-Item $page -ErrorAction SilentlyContinue
}
