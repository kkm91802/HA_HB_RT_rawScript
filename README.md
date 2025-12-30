# HA_HB_RT_rawScript

Vehicle harsh acceleration and braking detection using IMU MEMS (Bosch SMI230) data, with offline simulation based on recorded vehicle ECU CAN bus logs.

## 1. Overview

This project addresses the real-world problem of **driver behavior analysis**, specifically detecting **harsh acceleration and harsh braking events** using inertial sensor data.

The work was carried out during an **internship at Bosch Global Software Technologies Pvt. Ltd., Bangalore**, and focuses on replacing **lab-identified fixed threshold detection** with a **self-adaptive, signal-processing-based detection pipeline**.

The repository contains a **simulation-only, end-to-end pipeline** built on recorded vehicle data. It is **not a real-time ECU implementation**, but serves as a research and validation framework.

## 2. Motivation

Traditional driver behavior detection often relies on fixed thresholds calibrated in controlled lab environments. These thresholds:
- Do not generalize well across vehicles
- Are sensitive to sensor orientation and noise
- Fail under varying road and driving conditions

This project explores **adaptive and multi-method detection strategies** to improve robustness and reliability using real-world data.

## 3. Features & Capabilities

- End-to-end **offline simulation pipeline**
- IMU calibration using **three different calibration models**
- **Sensor orientation correction** for 6-axis IMU data
- Signal prediction using **Invariant Extended Kalman Filter (InEKF)**
- Signal processing techniques:
  - Butterworth low-pass filtering
  - Derivative-based detection
  - Energy envelope analysis
  - CUSUM
  - Continuous Wavelet Transform (CWT)
- **Adaptive thresholding** (self-adjusting detection logic)
- Polling / fusion logic to generate final detection output
- Python-based signal visualization and plotting
- Modular C-based implementation

Current detection focuses on **harsh acceleration and braking**. Rash turn detection is a planned extension.

## 4. Project Structure
├── src/ # Core C source files (including main.c)
├── include/ # Header files (.h)
├── scripts/ # Python scripts for plotting and analysis
├── csv/ # Processed IMU CSV data extracted from ECU logs
├── docs/ # Research papers, datasheets, references
├── CMakeLists.txt
└── README.md

## 5. Build Instructions

### Environment
- **Primary OS:** Ubuntu
- **Compiler:** GCC
- **Build system:** CMake
- **Dependencies:** Standard GCC / C libraries only

### Build Steps
mkdir build
cd build
cmake ..
make clean && make
./shared_folder

## 6. Execution Flow

The complete execution flow is implemented in main.c:

Convert recorded vehicle ECU logs into CSV format

Apply sensor orientation correction to align IMU data with vehicle frame

Perform IMU calibration

Predict calibrated signals using InEKF

Apply Butterworth low-pass filtering to remove high-frequency noise

Execute detection modules:

Adaptive thresholding

Derivative-based detection

Energy envelope

CUSUM

Continuous Wavelet Transform (CWT)

Poll detection outputs to generate final event decisions

Python scripts visualize and analyze the generated results

## 7. Data Handling

CSV files contain processed sensor data derived from vehicle ECU CAN logs

Original raw logs are confidential and not included

Data corresponds to Bosch SMI230 IMU

Provided CSVs are example datasets

Expected CSV Format

Users may replace CSV files using the following format:

time, ax, ay, az, gx, gy, gz

## 8. Cross-Platform Notes

Developed and validated primarily on Ubuntu

Windows was used only for:

Git version control

Repository management

.gitattributes is used for line-ending normalization (CRLF vs LF)

Windows runtime is not supported

Certain modules rely on Linux-specific functionality

Function replacements are required for Windows execution

## 9. Limitations & Future Work
Current Limitations

Offline simulation only

CSV-based file processing

Rash turn detection not fully implemented

Windows runtime unsupported

Limited dataset size

Future Work

Real-time ECU implementation

Removal of CSV parsing

Introduction of buffer-based streaming pipelines

Reduction/removal of library dependencies

Extension to gyro-based rash turn (left/right) detection

## 10. Internship & Disclaimer

This project was developed during an internship at:

Bosch Global Software Technologies Pvt. Ltd., Bangalore

Disclaimer

This repository represents independent academic/research work

It does not represent official Bosch software, IP, or products

The implementation is simulation-only and non-production

## 11. Author & Rights

Author:Kamal Kishore Majhi
GitHub: kkm91802

## Rights & License

All rights reserved.

No part of this repository may be used, modified, or redistributed without explicit permission from the author.
