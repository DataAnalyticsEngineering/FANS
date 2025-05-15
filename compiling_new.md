# 🔧 Compilation and Build Guide for FANS

This document outlines how to configure, build, and maintain reproducible builds of the **FANS** project using modern CMake presets, with support for multiple compilers and platforms.

---

## 🚀 Prerequisites

- **CMake ≥ 3.23** is required to support `CMakePresets.json` (version 6 or higher).
- Ensure a working C++ compiler is installed:
  - GCC (Linux)
  - Clang (Linux/macOS)
- CMake >= 3.23 can be installed from a package manager like `apt`, `brew`, or directly from [cmake.org](https://cmake.org/download/).

---

## ⚙️ Preset-Based Build System

CMake presets are defined in `CMakePresets.json` and allow standardized builds across configurations and platforms.

### ✅ Available Configure Presets

| Preset Name              | Compiler       | Platform | C++ Stdlib   | Notes                                                  |
|--------------------------|----------------|----------|--------------|--------------------------------------------------------|
| `gcc-release`            | GCC            | Linux    | `libstdc++`  | Standard Linux GCC release build                      |
| `gcc-debug`              | GCC            | Linux    | `libstdc++`  | Debug variant                                          |
| `clang-release-linux`    | Clang          | Linux    | `libstdc++`  | Requires GCC installed for `libstdc++` linking        |
| `clang-debug-linux`      | Clang          | Linux    | `libstdc++`  | Requires `libstdc++-dev` and `libomp` via `apt`       |
| `clang-release-macos`    | Clang (Apple)  | macOS    | `libc++`     | Uses macOS system Clang and default libc++            |
| `clang-debug-macos`      | Clang (Apple)  | macOS    | `libc++`     | Same as above with Debug mode                         |

### 🔨 Available Build Presets

| Preset Name              | Description                         |
|--------------------------|-------------------------------------|
| `build-gcc-release`      | Build using `gcc-release` preset    |
| `build-gcc-debug`        | Build using `gcc-debug` preset      |
| `build-clang-release-linux` | Build with Linux Clang release  |
| `build-clang-debug-linux`   | Build with Linux Clang debug    |
| `build-clang-release-macos` | Build with macOS Clang release  |
| `build-clang-debug-macos`   | Build with macOS Clang debug    |

---

## 🏗️ Building the Project

### 1. Clone the repository

```bash
git clone https://github.com/DataAnalyticsEngineering/FANS.git
cd FANS
```

### 2. Configure using CMake presets

To configure the build system, choose the appropriate preset for your platform and toolchain.

#### Example: GCC Release (Linux)

```bash
cmake --preset gcc-release
```

#### Example: GCC Debug (Linux)

```bash
cmake --preset gcc-debug
```

#### Example: Clang Release (Linux)

```bash
cmake --preset clang-release-linux
```

#### Example: Clang Debug (Linux)

```bash
cmake --preset clang-debug-linux
```

#### Example: Clang Release (macOS)

```bash
cmake --preset clang-release-macos
```

#### Example: Clang Debug (macOS)

```bash
cmake --preset clang-debug-macos
```

---

### 3. Build the Project

Use the following command to build the project after configuration:

```bash
cmake --build --preset <preset>
```

For example, if you used the `gcc-release` preset:

```bash
cmake --build --preset gcc-release
```

---

### Notes on Clang + libstdc++ (Linux)

If you're using Clang on Linux, and linking against `libstdc++`, you may run into issues where the linker cannot find `-lstdc++`. This is due to missing development libraries that Clang depends on for C++ standard support.
To make clang find those, there is a special `"CMAKE_CXX_FLAGS": "-stdlib=libstdc++"` in the CMAKE-Preset.

Ensure you have installed the required GCC development libraries:

```bash
sudo apt install g++ libstdc++-<version>-dev
```

Where `<version>` should match your system's default GCC version. You can discover available versions via:

```bash
apt search libstdc++
```

---

### Platform Distinctions in Presets

Because Clang on Linux often requires explicit linking to `libstdc++`, the presets distinguish between:

- `clang-release-linux`, `clang-debug-linux` – for Linux builds using Clang with system-wide libstdc++.
- `clang-release-macos`, `clang-debug-macos` – for macOS builds using Clang with system libc++.

This allows for reproducible builds on multiple platforms with appropriate linking and standard library settings.
