![screenshot](res/chessboard.png)

# `lensprofiler`

`lensprofiler` application implements a workflow to create lens profiles, that can be used together with [MRTech IFF SDK](https://mr-te.ch/iff-sdk).
It shares most of the C++ code with [`imagebroker`](https://github.com/mr-technologies/imagebroker) example IFF SDK application, but also includes `ldc.py` Python script for chessboard target auto-detection and radial lens distortion model estimation, which is based on the source code provided with [Algebraic Lens Distortion Model Estimation](https://www.ipol.im/pub/art/2011/ags-alde/) article by Luis Alvarez, Luis Gomez, and J. Rafael Sendra published in [Image Processing On Line](https://www.ipol.im/) journal.
Application is located in `samples/06_lens` directory of IFF SDK package.
It comes with example configuration files (`lensprofiler.json` and `res/ldc.json`) suited for XIMEA cameras and standard 10x7 chessboard target (such as included `res/chessboard.png`).
See `linux` and `windows` directories for helper scripts to install required dependencies (e.g. [OpenCV](https://opencv.org/) library).
Operation is controlled using a keyboard:

* `1` decreases exposure
* `2` increases exposure
* `Tab` captures an image and starts the profile generation procedure (further instructions are shown on the screen)
