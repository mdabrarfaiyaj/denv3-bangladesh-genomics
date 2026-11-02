#!/usr/bin/env bash
#
# 14_system_check.sh
# ------------------
# Can this machine actually run a containerised Nextflow pipeline?
#
# Checks facts only and prints them plainly. It installs nothing and changes
# nothing. The point is to find out BEFORE committing to an approach whether
# the approach fits the machine - 5 GB of RAM and Docker Desktop are not an
# obvious combination.
#
# What matters and why:
#   WSL version   Docker needs WSL2. WSL1 cannot run it.
#   RAM           Nextflow's JVM, Docker's VM and IQ-TREE all want memory at
#                 the same time. The whole-genome tree alone used most of a
#                 4-core machine for ten minutes.
#   Java          Nextflow requires Java 17-24 (21 is the safe choice).
#   Docker        Either Docker Desktop with WSL integration, or the Docker
#                 engine installed inside WSL2 directly (needs systemd).
#   Filesystem    Work under ~ (ext4), never under /mnt/c or /mnt/f. Windows
#                 drives through the 9p bridge are 10-50x slower, which turns a
#                 ten-minute pipeline into an hour.
#   Network       Container images come from quay.io and docker.io.
#
# Run from the project folder:  bash scripts/14_system_check.sh

echo "=============================================================="
echo "  SYSTEM CHECK - facts only, nothing is installed or changed"
echo "=============================================================="

say()  { printf "  %-26s %s\n" "$1" "$2"; }
have() { command -v "$1" >/dev/null 2>&1; }

# ------------------------------------------------------------------ 1. OS --
echo
echo "-- 1. Operating system --"
say "kernel" "$(uname -r)"
if grep -qi "microsoft" /proc/version 2>/dev/null; then
  if uname -r | grep -qi "WSL2\|microsoft-standard-WSL2"; then
    say "WSL version" "WSL2  (Docker is possible)"
  else
    say "WSL version" "WSL1 or unclear  -> Docker will NOT work; see verdict"
  fi
else
  say "WSL" "not WSL (native Linux)"
fi
[ -f /etc/os-release ] && say "distro" "$(. /etc/os-release; echo "$PRETTY_NAME")"
if have systemctl && systemctl is-system-running >/dev/null 2>&1; then
  say "systemd" "running  (Docker engine can run inside WSL)"
else
  say "systemd" "not running  (needed only for Docker engine inside WSL)"
fi

# --------------------------------------------------------------- 2. CPU/RAM --
echo
echo "-- 2. CPU and memory --"
say "CPU cores" "$(nproc 2>/dev/null || echo '?')"
have lscpu && say "CPU model" "$(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')"
if [ -r /proc/meminfo ]; then
  memkb=$(awk '/MemTotal/{print $2}' /proc/meminfo)
  swapkb=$(awk '/SwapTotal/{print $2}' /proc/meminfo)
  memgb=$(awk -v k="$memkb" 'BEGIN{printf "%.1f", k/1048576}')
  say "total RAM" "${memgb} GB"
  say "total swap" "$(awk -v k="$swapkb" 'BEGIN{printf "%.1f", k/1048576}') GB"
  say "available now" "$(awk '/MemAvailable/{printf "%.1f GB", $2/1048576}' /proc/meminfo)"
fi
# Windows-side limit, if the user set one
for u in /mnt/c/Users/*/; do
  [ -f "${u}.wslconfig" ] && { say ".wslconfig found" "${u}.wslconfig"; grep -E "memory|processors|swap" "${u}.wslconfig" 2>/dev/null | sed 's/^/                             /'; }
done

# ------------------------------------------------------------------ 3. Disk --
echo
echo "-- 3. Disk and filesystem --"
say "project path" "$(pwd)"
case "$(pwd)" in
  /mnt/*) say "filesystem" "WINDOWS DRIVE (9p) - SLOW, move the project under ~" ;;
  *)      say "filesystem" "$(df -T . 2>/dev/null | awk 'NR==2{print $2}') (native Linux - good)" ;;
esac
df -h . 2>/dev/null | awk 'NR==2{printf "  %-26s %s free of %s\n", "space here", $4, $2}'
[ -d results ] && say "results/ size" "$(du -sh results 2>/dev/null | cut -f1)"
[ -d data ] && say "data/ size" "$(du -sh data 2>/dev/null | cut -f1)"

# ------------------------------------------------------------------ 4. Java --
echo
echo "-- 4. Java (Nextflow needs 17-24) --"
if have java; then
  # java prints JAVA_TOOL_OPTIONS and other noise before the version, so match
  # the version line itself rather than taking the first line
  v=$(java -version 2>&1 | grep -m1 'version "')
  [ -z "$v" ] && v=$(java -version 2>&1 | head -1)
  say "java" "$v"
  major=$(printf '%s' "$v" | sed -E 's/.*version "([0-9]+).*/\1/')
  # "1.8.0_x" style means Java 8
  [ "$major" = "1" ] && major=8
  if [ -n "$major" ] && [ "$major" -ge 17 ] 2>/dev/null && [ "$major" -le 24 ] 2>/dev/null; then
    say "verdict" "OK for Nextflow"
  else
    say "verdict" "version may not suit Nextflow (wants 17-24)"
  fi
else
  say "java" "NOT INSTALLED  -> Nextflow cannot run"
fi

# -------------------------------------------------------------- 5. Nextflow --
echo
echo "-- 5. Nextflow --"
if have nextflow; then
  say "nextflow" "$(nextflow -v 2>&1 | head -1)"
else
  say "nextflow" "not on PATH (it is in environment.yml; conda activate denv3)"
fi

# --------------------------------------------------------- 6. Containers ----
echo
echo "-- 6. Containers --"
if have docker; then
  say "docker client" "$(docker --version 2>&1)"
  if docker info >/dev/null 2>&1; then
    say "docker daemon" "RUNNING"
    say "docker root" "$(docker info --format '{{.DockerRootDir}}' 2>/dev/null)"
    say "docker memory" "$(docker info --format '{{.MemTotal}}' 2>/dev/null | awk '{printf "%.1f GB", $1/1073741824}')"
  else
    say "docker daemon" "NOT reachable (Docker Desktop not started, or no WSL integration)"
  fi
else
  say "docker" "not installed"
fi
for s in singularity apptainer; do
  have $s && say "$s" "$($s --version 2>&1 | head -1)"
done
if ! have docker && ! have singularity && ! have apptainer; then
  say "verdict" "no container runtime - conda profile is the fallback"
fi

# ------------------------------------------------------------- 7. Toolchain --
echo
echo "-- 7. Other tools --"
for t in conda mamba git gh python Rscript iqtree iqtree2 treetime nextclade mafft; do
  have $t && say "$t" "$(command -v $t)"
done

# ---------------------------------------------------------------- 8. Network --
echo
echo "-- 8. Network reachability (container registries and code hosts) --"
for host in quay.io registry-1.docker.io github.com anaconda.org; do
  if curl -s -o /dev/null -m 8 -w "%{http_code}" "https://$host" 2>/dev/null | grep -qE "^(2|3|4)"; then
    say "$host" "reachable"
  else
    say "$host" "NOT reachable or timed out"
  fi
done

# ----------------------------------------------------------------- verdict --
echo
echo "=============================================================="
echo "  WHAT THIS MEANS"
echo "=============================================================="
memgb_int=$(awk '/MemTotal/{printf "%d", $2/1048576}' /proc/meminfo 2>/dev/null || echo 0)
cores=$(nproc 2>/dev/null || echo 0)

echo
echo "  cores=$cores  RAM=${memgb_int}GB"
if [ "$memgb_int" -lt 8 ]; then
cat <<'EOF'

  RAM is the binding constraint. Under 8 GB, running Nextflow (a JVM),
  Docker (a second Linux VM under Windows) and IQ-TREE at the same time
  will be tight. That does NOT stop you building a proper pipeline - it
  changes which profile you run locally:

      local runs  -> -profile conda   (no container overhead)
      CI runs     -> -profile docker  (GitHub Actions gives you a clean
                                       runner with Docker preinstalled)

  The repository still ships a docker profile and is still containerised
  and reproducible, which is what matters for the portfolio. It simply is
  not the profile you use on this laptop.
EOF
else
cat <<'EOF'

  Enough memory to run the docker profile locally as well as in CI.
EOF
fi
cat <<'EOF'

  Things to fix before starting, if the checks above flagged them:
    * project NOT under ~ (e.g. under /mnt/c or /mnt/f) - move it. Windows
      drives via the 9p bridge are many times slower for the thousands of
      small file operations a pipeline makes.
    * no Java - install it (e.g. conda install -c conda-forge openjdk=21).
    * Docker daemon unreachable - either start Docker Desktop and enable WSL
      integration for this distro, or skip Docker locally and use conda.

  Nothing here blocks the work. The point of checking first is to pick the
  profile that fits the machine rather than discovering it halfway through.
EOF
