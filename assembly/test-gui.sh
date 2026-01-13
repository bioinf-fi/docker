#!/bin/bash
# Test script to verify X11/GUI setup for Docker
# Run this BEFORE starting the container to diagnose issues

set -e

echo "=========================================="
echo "X11/GUI Setup Verification"
echo "=========================================="
echo ""

# Test 1: Check if DISPLAY is set
echo "Test 1: Checking DISPLAY environment variable..."
if [ -z "$DISPLAY" ]; then
    echo "  ❌ FAIL: DISPLAY is not set"
    echo "  → You may not be on a graphical desktop"
    exit 1
else
    echo "  ✓ PASS: DISPLAY=$DISPLAY"
fi
echo ""

# Test 2: Check X11 socket
echo "Test 2: Checking X11 socket..."
if [ -S "/tmp/.X11-unix/X${DISPLAY#*:}" ] || [ -S "/tmp/.X11-unix/X0" ]; then
    echo "  ✓ PASS: X11 socket found"
else
    echo "  ⚠ WARNING: X11 socket not found at expected location"
    echo "  → This is normal on Wayland systems"
fi
echo ""

# Test 3: Check .Xauthority
echo "Test 3: Checking .Xauthority file..."
if [ -f "$HOME/.Xauthority" ]; then
    echo "  ✓ PASS: .Xauthority exists at $HOME/.Xauthority"
    XAUTH_ENTRIES=$(xauth list 2>/dev/null | wc -l || echo "0")
    if [ "$XAUTH_ENTRIES" -gt 0 ]; then
        echo "  ✓ PASS: .Xauthority contains $XAUTH_ENTRIES entries"
    else
        echo "  ⚠ WARNING: .Xauthority is empty or unreadable"
    fi
else
    echo "  ⚠ WARNING: .Xauthority not found"
    echo "  → Will use relaxed X11 access mode"
fi
echo ""

# Test 4: Check if xhost is available
echo "Test 4: Checking xhost utility..."
if command -v xhost >/dev/null 2>&1; then
    echo "  ✓ PASS: xhost is available"
    # Try to check current access control
    if xhost 2>&1 | grep -q "access control enabled"; then
        echo "  ℹ INFO: X11 access control is enabled (secure)"
    elif xhost 2>&1 | grep -q "access control disabled"; then
        echo "  ⚠ WARNING: X11 access control is disabled (insecure but will work)"
    fi
else
    echo "  ⚠ WARNING: xhost not found"
    echo "  → You may need to install xhost for X11 access control"
fi
echo ""

# Test 5: Check container runtime (Docker or Podman)
echo "Test 5: Checking container runtime..."
CONTAINER_RUNTIME=""
if command -v podman >/dev/null 2>&1; then
    CONTAINER_RUNTIME="podman"
    echo "  ✓ PASS: Podman is installed"
    if podman ps >/dev/null 2>&1; then
        echo "  ✓ PASS: Podman is running and accessible"
    else
        echo "  ❌ FAIL: Podman is not accessible"
        exit 1
    fi
elif command -v docker >/dev/null 2>&1; then
    CONTAINER_RUNTIME="docker"
    echo "  ✓ PASS: Docker is installed"
    if docker ps >/dev/null 2>&1; then
        echo "  ✓ PASS: Docker is running and accessible"
    else
        echo "  ❌ FAIL: Docker is not running or you don't have permission"
        echo "  → Try: sudo usermod -aG docker $USER (then logout/login)"
        exit 1
    fi
else
    echo "  ❌ FAIL: Neither Docker nor Podman is installed"
    echo "  → Please install Docker or Podman"
    exit 1
fi
echo "  ℹ INFO: Using $CONTAINER_RUNTIME as container runtime"
echo ""

# Test 6: Try a simple X11 test
echo "Test 6: Testing X11 display connection..."
if command -v xdpyinfo >/dev/null 2>&1; then
    if xdpyinfo >/dev/null 2>&1; then
        echo "  ✓ PASS: X11 display is accessible"
        DISPLAY_DIM=$(xdpyinfo | grep dimensions | awk '{print $2}')
        echo "  ℹ INFO: Display dimensions: $DISPLAY_DIM"
    else
        echo "  ❌ FAIL: Cannot connect to X11 display"
        exit 1
    fi
else
    echo "  ⚠ WARNING: xdpyinfo not available (install x11-utils to test)"
fi
echo ""

# Test 7: Check display server type
echo "Test 7: Detecting display server..."
if [ "$XDG_SESSION_TYPE" = "wayland" ]; then
    echo "  ℹ INFO: Running on Wayland (with XWayland)"
    echo "  → GUI apps should work via XWayland compatibility"
elif [ "$XDG_SESSION_TYPE" = "x11" ]; then
    echo "  ✓ PASS: Running on native X11"
else
    echo "  ℹ INFO: Display server type unknown (likely X11)"
fi
echo ""

# Summary
echo "=========================================="
echo "Summary"
echo "=========================================="
echo ""
echo "Your system appears ready for containerized GUI applications!"
echo "Container runtime: $CONTAINER_RUNTIME"
echo ""
echo "Next steps:"
echo "  1. Build the image: make build"
echo "  2. Start with GUI:   make gui"
echo "  3. If step 2 fails:  make gui-simple"
echo ""
echo "Inside the container, run:"
echo "  - igv       (Integrative Genomics Viewer)"
echo "  - Bandage   (Assembly graph viewer)"
echo ""
