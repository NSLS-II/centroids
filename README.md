# centroids

## Build Status


| Build Type | Result |
|---|:-----:|
| **Python Linux Builds** | [![Build Status](https://dev.azure.com/nsls-ii/centroids/_apis/build/status/NSLS-II.centroids?branchName=master&jobName=BuildPython)](https://dev.azure.com/nsls-ii/centroids/_build/latest?definitionId=3&branchName=master) |
| **CXX Lib Linux Builds** | [![Build Status](https://dev.azure.com/nsls-ii/centroids/_apis/build/status/NSLS-II.centroids?branchName=master&jobName=BuildLib)](https://dev.azure.com/nsls-ii/centroids/_build/latest?definitionId=3&branchName=master) |
| **Documentation Build** | [![Build Status](https://dev.azure.com/nsls-ii/centroids/_apis/build/status/NSLS-II.centroids?branchName=master&jobName=BuildDocs)](https://dev.azure.com/nsls-ii/centroids/_build/latest?definitionId=3&branchName=master) |

## Introduction

Routines for performing single photon counting and centroiding of CCD data

## Documentation

The documentation can be found [online](https://nsls-ii.github.io/centroids/).

## Changes

### v0.2.0rc1

Release candidate for v0.2.

* Changed `MANIFEST.in` to allow for building with `pip`.

### v0.1.11

* Fixed bug in standard error for 2D fits calculations. Incorrect number of data points was used.
* Added range option to fitting parameters for _sigma_ and _position_
* Modified CI to report tests and coverage
* Improved unit tests to cover parameter range selection
* Added debug output to unit tests
* Updated docs