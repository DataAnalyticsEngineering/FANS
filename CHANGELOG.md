# FANS Changelog

## latest

- Relaxed versions in Pixi build workflow [#86](https://github.com/DataAnalyticsEngineering/FANS/pull/86)
- Bump version in MacOS CI [#85](https://github.com/DataAnalyticsEngineering/FANS/pull/85)
- Pixi build workflow for MacOS and Linux [#82](https://github.com/DataAnalyticsEngineering/FANS/pull/82) [#84](https://github.com/DataAnalyticsEngineering/FANS/pull/84)
- Introduce conftest.py for pytest-fixture [#82](https://github.com/DataAnalyticsEngineering/FANS/pull/82)
- Use dashboard environment for regular pytest [#81](https://github.com/DataAnalyticsEngineering/FANS/pull/81)
- Fix bug in J2 plasticity routine [#80](https://github.com/DataAnalyticsEngineering/FANS/pull/80)

## v0.4.2

- Reduce dependencies of dashboard for cf-package https://github.com/DataAnalyticsEngineering/FANS/pull/71
- Add pixi task `h52xdmf` to generate XDMF from H5 files directly as `pixi run h52xdmf {h5filepath}`
- Add logo of FANS to the README https://github.com/DataAnalyticsEngineering/FANS/pull/70

## v0.4.1

- remove std::sqrt from constexpr - failed on Clang https://github.com/DataAnalyticsEngineering/FANS/pull/64

## v0.4.0

- Support compilaion on MacOS X via conda-forge https://github.com/DataAnalyticsEngineering/FANS/pull/59
- Add support for macroscale mixed stress-strain boundary conditions https://github.com/DataAnalyticsEngineering/FANS/pull/58
- Add grain boundary diffusion material model for polycrystals https://github.com/DataAnalyticsEngineering/FANS/pull/52
- Add a pixi environment for the FANS_dashboard and some tests https://github.com/DataAnalyticsEngineering/FANS/pull/55
- Remove MPI initialization from pyFANS and add an integration test for it https://github.com/DataAnalyticsEngineering/FANS/pull/46
- Native support for MacOS https://github.com/DataAnalyticsEngineering/FANS/pull/25
- Remove Ubuntu 20.04 from testing and Docker support https://github.com/DataAnalyticsEngineering/FANS/pull/51
- Add support for `--version` command line argument for checking the version of FANS
- Modify way to provide micro structure in JSON input https://github.com/DataAnalyticsEngineering/FANS/pull/43
- Add conda package for FANS https://github.com/DataAnalyticsEngineering/FANS/pull/39
- Introduce system for checking compiler flags: `avx2` & `fma` https://github.com/DataAnalyticsEngineering/FANS/pull/34
- Add `results_prefix` field in the JSON input https://github.com/DataAnalyticsEngineering/FANS/pull/36
- Build FANS as a library to be coupled to a macro-scale simulation via preCICE and the Micro Manager https://github.com/DataAnalyticsEngineering/FANS/pull/23

## v0.3.0

- Added Linear thermal and mechanical triclinic material models https://github.com/DataAnalyticsEngineering/FANS/pull/32
- Added API to get homogenized stress and homogenized tangent https://github.com/DataAnalyticsEngineering/FANS/pull/31

## v0.2.0

- Add integration tests https://github.com/DataAnalyticsEngineering/FANS/pull/20
- Add GitHub Action workflow to build and test FANS https://github.com/DataAnalyticsEngineering/FANS/pull/19

## v0.1.2

- Update TIK GitHub links in the documentation to public GitHub links https://github.com/DataAnalyticsEngineering/FANS/pull/13

## v0.1.1

- Disable sorting of includes in clang-format https://github.com/DataAnalyticsEngineering/FANS/pull/7

## v0.1.0

- Add release guide and a Changelog file https://github.com/DataAnalyticsEngineering/FANS/pull/4
- Add clang-format check and format all relevant files https://github.com/DataAnalyticsEngineering/FANS/pull/1
