#!/usr/bin/env bash
# Install (or refresh) the weekly GitHub traffic collector as a launchd agent.
#
# macOS denies a launchd agent access to ~/Desktop, so the collector runs from
# ~/Library/Application Support/smart-traffic and the repository's metrics/ is a symlink
# to the data there. Re-run this after changing utils/collect_github_traffic.py.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
APP="$HOME/Library/Application Support/smart-traffic"
AGENT="$HOME/Library/LaunchAgents/com.monkiky.smart-traffic.plist"
LABEL="com.monkiky.smart-traffic"

mkdir -p "$APP/metrics" "$HOME/Library/LaunchAgents"
cp -f "$REPO/utils/collect_github_traffic.py" "$APP/"
# the plist in the repository is a template: no machine-specific paths are committed
sed -e "s|__APP_DIR__|$APP|g" -e "s|__REPO__|$REPO|g" \
    -e "s|/opt/homebrew/bin|$(dirname "$(command -v gh)"):/opt/homebrew/bin|" \
    "$REPO/utils/launchd/$LABEL.plist" > "$AGENT"

# the repository sees the data through a symlink, so a manual run and the scheduled run
# write to the same files
if [[ ! -L "$REPO/metrics" ]]; then
    [[ -d "$REPO/metrics" ]] && cp -pn "$REPO/metrics"/*.csv "$APP/metrics/" 2>/dev/null || true
    rm -rf "$REPO/metrics"
    ln -s "$APP/metrics" "$REPO/metrics"
fi

launchctl bootout "gui/$(id -u)/$LABEL" 2>/dev/null || true
launchctl bootstrap "gui/$(id -u)" "$AGENT"
echo "installed: $LABEL (Mondays 10:30)"
echo "  run now:  launchctl kickstart -p gui/$(id -u)/$LABEL"
echo "  log:      $APP/collector.log"
echo "  remove:   launchctl bootout gui/$(id -u)/$LABEL && rm '$AGENT'"
