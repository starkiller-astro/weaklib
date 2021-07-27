# Changelog
All notable changes to this project will be documented in this file.

## [Unreleased]

## [1.0.1] - 2021-07-23
### Added
- Added 'wlInterpolationUtilitiesModule.F90' for basic interpolation routines, such as 'Index1D', 'LinearInterp_Array_Point'.
- Added 'wlOpacityPerformanceTest.F90' as an unit test.

### Changed
- Cleaned up `wlInterpolationModule.F90` and moved utility functions. Added OpenMP directives that were missing.

### Removed
- Removed OpacityAlignTest and OpacityTableResolutionTest from standard unit test list.

## [1.0.0] - 2020-07-?
### Added
- New README 
- New "Guide for Developers"
- New Changelog

### Changed
- Moved 'ComputeTempFromX' routines from Library/wlInterpolationModule.F90 to EOSSource/wlEOSInterpolationModule.f90
- Cosmetic changes in Library/wlInterpolationModule.F90

### Removded
