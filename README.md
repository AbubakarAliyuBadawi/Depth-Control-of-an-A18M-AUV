# Depth Control of A18M Autonomous Underwater Vehicle

## Project Overview
This repository contains the MATLAB and Simulink models used to develop and test a Linear Quadratic Regulator (LQR) for controlling the A18M Autonomous Underwater Vehicle (AUV) depth. The project's goal was to optimize depth control with minimal pitch oscillations, crucial for the AUV's mission of underwater mine detection

## Group Members
- **[Abubakar Aliyu Badawi](https://www.linkedin.com/in/abubakar-aliyu-badawi-41533b150/)**
- **Krittapat Onthuam**
- **Tayyab Tahir**
- **Ishfaq Bhat**

<div align="center">
  <img src="Plots/A-18%20AUV.png" alt="AUV Model" width="400" height="300">
  <p><em>Figure 1: A-18 Autonomous Underwater Vehicle used in simulations.</em></p>
</div>

### Background
The A18M, an AUV developed by ECARobotics, is designed for efficient and accurate underwater mine detection. Effective mine detection requires a stable platform, necessitating precise control over the vehicle's movements, especially in the vertical plane.

### Objectives
- **Replace Traditional PID Controllers**: Implement an LQR controller for improved response and stability.
- **Enhance Depth Control**: Achieve precise depth targeting with reduced oscillations and improved response times.
- **Simulation and Validation**: Use MATLAB and Simulink to simulate the controller's performance and validate its effectiveness.

## Installation

### Prerequisites
- MATLAB (2020 or later recommended)
- Simulink
- MATLAB Control System Toolbox

### Setup
Clone this repository to your local machine using:
```bash
git clone http://github.com/yourusername/auv-depth-control.git
cd auv-depth-control

