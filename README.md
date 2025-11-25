# Enhanced in vivo thermoacoustic imaging of intracerebral hemorrhage in rat models through an optimal microwave illumination strategy

This repository contains the MATLAB code used for Thermoacoustic Tomography (TAT) image reconstruction and enhancement in our study on intracerebral hemorrhage (ICH) imaging. The implementation includes the Delay-and-Sum (DAS) algorithm for core image formation and the Multi-Scale Retinex-baesd algorithm with logarithmic depth weighting (Log-MSR). for image quality improvement. The code is designed to process raw data obtained from both normal and ICH rat models, with and without the metal reflector illumination strategy.

**Data & Workflow:**
- **Datasets:** Four `.lvm` files are provided, comprising normal and lesion models, each with and without the metal reflector, the number in parentheses represents the image reconstruction radius, R0.
- **Reconstruction:** Reconstruct the image by running the `DAS.m` code in MATLAB and save the output as a `.fig` file.
- **Enhancement:** Run `Log+MSR.m` and select the saved `.fig` file to apply the image enhancement.

This work is associated with the paper **"Enhanced in vivo thermoacoustic imaging of intracerebral hemorrhage in rat models through an optimal microwave illumination strategy"** conducted at the **Chongqing University of Posts and Telecommunications (CQUPT)**.

THIS CODE IS PROVIDED FOR ACADEMIC RESEARCH PURPOSES ONLY, ON AN "AS IS" BASIS. THE AUTHOR MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING ITS COMPLETENESS, ACCURACY, OR FITNESS FOR A PARTICULAR PURPOSE. USERS ASSUME ALL RISKS ASSOCIATED WITH ITS USE. THE AUTHOR SHALL NOT BE LIABLE FOR ANY LOSS OR DAMAGE ARISING FROM THE USE OF THIS CODE.IF YOU ARE PUBLISHING ANY WORK WHERE THIS CODE OR A VARIANT OF IT HAS BEEN USED, PLEASE REMEMBER THAT IT WAS OBTAINED FREE OF CHARGE. WE KINDLY ASK THAT YOU REFERENCE OUR PAPER IN THE PUBLICATION. 
