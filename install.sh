#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# ==============================================================================
# Biologist-Friendly One-Click Installer
# Installs: minimap2, r2rtr, raEDClust, pyfastx, and k8 (for srfutils.js)
# ==============================================================================

# --- 0. Configuration & Environment Checks ---

# Define Colors for Output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Starting Installation Pipeline ===${NC}"

# Get current directory
ROOT_DIR=$(pwd)

# Detect CPU cores for parallel compilation (make -j)
if command -v nproc > /dev/null; then
    THREADS=$(nproc)
elif [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    THREADS=$(sysctl -n hw.ncpu)
else
    THREADS=4
fi
echo -e "Detected ${THREADS} CPU cores. Will use ${THREADS} threads for compilation."

# Check basic system dependencies
echo -e "${YELLOW}[Check] Verifying system dependencies...${NC}"
for cmd in make gcc g++ python3 pip3 curl tar; do
    if ! command -v $cmd &> /dev/null; then
        echo -e "${RED}Error: Command '$cmd' not found. Please install it first.${NC}"
        exit 1
    fi
done

# --- 1. Install Minimap2 ---
echo -e "${YELLOW}[Step 1/5] Installing Minimap2...${NC}"
if [ -d "minimap2" ]; then
    cd minimap2
    make -j"$THREADS"
    cd "$ROOT_DIR"
    echo -e "${GREEN}✔ Minimap2 installed successfully.${NC}"
else
    echo -e "${RED}Error: Directory 'minimap2' not found! Please ensure you are in the project root.${NC}"
    exit 1
fi

# --- 2. Install r2rtr ---
echo -e "${YELLOW}[Step 2/5] Installing r2rtr...${NC}"
if [ -d "r2rtr" ]; then
    cd r2rtr
    make -j"$THREADS"
    cd "$ROOT_DIR"
    echo -e "${GREEN}✔ r2rtr installed successfully.${NC}"
else
    echo -e "${RED}Error: Directory 'r2rtr' not found!${NC}"
    exit 1
fi

# --- 3. Install raEDClust ---
echo -e "${YELLOW}[Step 3/5] Installing raEDClust...${NC}"
if [ -d "raEDClust" ]; then
    cd raEDClust
    make -j"$THREADS"
    cd "$ROOT_DIR"
    echo -e "${GREEN}✔ raEDClust installed successfully.${NC}"
else
    echo -e "${RED}Error: Directory 'raEDClust' not found!${NC}"
    exit 1
fi

# --- 4. Install Python Dependencies (pyfastx) ---
echo -e "${YELLOW}[Step 4/5] Installing Python dependencies (pyfastx)...${NC}"
# Using --user to avoid permission issues without sudo
if pip3 install pyfastx --user; then
    echo -e "${GREEN}✔ pyfastx installed successfully.${NC}"
else
    echo -e "${RED}Error: Failed to install pyfastx. Check your internet connection or pip settings.${NC}"
    exit 1
fi

# --- 5. Install k8 (Javascript Shell) ---
echo -e "${YELLOW}[Step 5/5] Installing k8 (Javascript shell)...${NC}"

# Prepare installation path $HOME/bin (Standard user binary path)
INSTALL_BIN="$HOME/bin"
mkdir -p "$INSTALL_BIN"

# Define version and URL
K8_VER="0.2.4"
K8_URL="https://github.com/attractivechaos/k8/releases/download/v${K8_VER}/k8-${K8_VER}.tar.bz2"

echo "Downloading k8 from $K8_URL ..."
if curl -L "$K8_URL" | tar -jxf - ; then
    # Detect OS (Linux vs Mac)
    OS_NAME=$(uname -s)
    if [[ "$OS_NAME" == "Darwin" ]]; then
        K8_BINARY="k8-${K8_VER}/k8-Darwin"
    else
        K8_BINARY="k8-${K8_VER}/k8-Linux"
    fi

    # Copy binary to $HOME/bin
    if [ -f "$K8_BINARY" ]; then
        cp "$K8_BINARY" "$INSTALL_BIN/k8"
        chmod +x "$INSTALL_BIN/k8"
        rm -rf "k8-${K8_VER}" # Clean up temporary files
        echo -e "${GREEN}✔ k8 installed to $INSTALL_BIN/k8${NC}"
    else
        echo -e "${RED}Error: Could not find k8 binary for $OS_NAME inside the downloaded package.${NC}"
        exit 1
    fi
else
    echo -e "${RED}Error: Failed to download k8. Please check your internet connection.${NC}"
    exit 1
fi

# --- Completion & Environment Setup ---
echo -e "\n${GREEN}=== Installation Completed Successfully! ===${NC}"
echo ""
echo -e "Please ensure ${YELLOW}$HOME/bin${NC} is in your PATH to run srfutils.js."
echo -e "You can add the following line to your ~/.bashrc or ~/.zshrc:"
echo ""
echo -e "    ${GREEN}export PATH=\"\$HOME/bin:\$PATH\"${NC}"
echo ""
echo "After adding, run 'source ~/.bashrc' or restart your terminal."