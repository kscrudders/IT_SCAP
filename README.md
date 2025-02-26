<a id="readme-top"></a>
<!--
*** Thanks for checking out the IT_SCAP-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** 
*** I imagine a world where scientific knowledge provides solutions for every health challenge, enabling everyone to live with autonomy, freedom, and well-being.
*** I created this project so that I might streamline taking raw microscopy data in my PhD and convert that in biological insights that might aid understanding the next generation of engineered T cell immunotherapies.
*** I hope this could be useful to a few future scienctist in whatever pursuit they are taking on. 
*** I would be overjoyed to help enable you to make discoveries and share knowlegde with humanity.
-->

<!-- PROJECT LOGO --> <br /> <div align="center">   <a href="https://github.com/kscrudders/IT_SCAP"> <img src="images/IT_SCAP_projectlogo.png" alt="Logo" width="300" height="100"> </a> <h3 align="center">IRM/TIRF Single-Cell semiAutomated Pipeline (IT_SCAP)</h3> <p align="center"> A MATLAB-based workflow for segmenting and analyzing cell ROIs in TIRF/IRM datasets <br /> <a href="https://github.com/your_username/IT_SCAP"><strong>Explore the docs »</strong></a> <br /> <br /> <a href="https://github.com/kscrudders/IT_SCAP/issues">Report Bug</a> · <a href="https://github.com/kscrudders/IT_SCAP/issues">Request Feature</a> </p> </div> <!-- TABLE OF CONTENTS --> <details> <summary>Table of Contents</summary> <ol> <li><a href="#about-the-project">About The Project</a></li> <li><a href="#built-with">Built With</a></li> <li><a href="#getting-started">Getting Started</a> <ul> <li><a href="#prerequisites">Prerequisites</a></li> <li><a href="#installation">Installation</a></li> </ul> </li> <li><a href="#usage">Usage</a></li> <li><a href="#roadmap">Roadmap</a></li> <li><a href="#contributing">Contributing</a></li> <li><a href="#license">License</a></li> <li><a href="#contact">Contact</a></li> <li><a href="#acknowledgments">Acknowledgments</a></li> </ol> </details> <!-- ABOUT THE PROJECT -->
About The Project
IT_SCAP organizes and streamlines importing, processing, segmenting, and analyzing IRM paired TIRF microscopy data. It helps users:

Organize image data and metadata for multiple channels.
Perform background, shade, and bleaching corrections.
Define regions of interest (ROIs) for single cells.
Save processed data and ROI-based sub-images for further analysis.
This script is particularly useful as a pipiline to anyone performing single-cell time acquisitions with IRM imaging to segment cells and TIRF imaging for biological readouts.
This script is currently optimized for use with Nikon microscopes (.nd2 image files). 
TIF image files are also acceptable file inputs (format as interleaved channels in time. [i.e. ch1 time1, ch2 time1, ch1 time 2, ch2 time 2, etc])

<p align="right">(<a href="#readme-top">back to top</a>)</p> <!-- BUILT WITH -->

Built With Matlab and Matlab Toolboxs:
* Matlab vR2024b
* Computer Vision Toolbox
* Curve Fitting Toolbox
* Image Processing Toolbox
* Statistics and Machine Learning Toolbox
	
External Supporting Matlab scripts:
* bfopen - https://www.openmicroscopy.org/bio-formats/downloads/
* getkey - https://www.mathworks.com/matlabcentral/fileexchange/7465-getkey
* createMask, rescale, rdouble - https://github.com/DanuserLab/u-track
* ND2info - https://www.mathworks.com/matlabcentral/fileexchange/73271-nd2-reader
* plotSpread - https://www.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot

<p align="right">(<a href="#readme-top">back to top</a>)</p> <!-- GETTING STARTED -->

Getting Started </p>
These instructions will guide you in preparing your local environment to run IT_SCAP.m.

Prerequisites:
* MATLAB (R2019b or newer recommended)
* Above toolboxs
* Above external matlab scripts 
* Shade Image data
* Gain calibration data for EMCCD camera

Installation:
* Clone the repository:
* sh
* Copy
* Edit
* git clone https://github.com/kscrudders/IT_SCAP.git

Add to MATLAB path:
* Open MATLAB.
* Go to Home > Set Path > Add with Subfolders and select the cloned folder.

<p align="right">(<a href="#readme-top">back to top</a>)</p> <!-- USAGE EXAMPLES -->

Usage:
* In MATLAB, open IT_SCAP.m.
* Set file paths under Section_01 (Import Data). Adjust data_dir and save_data_in_this_folder to your local directories where raw image data and output will be stored.
Adjust acquisition parameters:
* Define channel imaging frequencies, and select gain, shade, and/or bleach correction setting, and so on in Sections_02–03.
Select and edit ROIs, to isolate cells.
* Remaining sections should be executed one by one, so you can inspect intermediate outputs and fix any issues (e.g., ROI adjustments, segmentation, tracking, etc).
Export: The script can produce TIF files, and .mat files containing processed data and ROI-based extractions.
* For step-by-step details, see the comments within the script. And/or walk through example on YouTube: [WIP]

<p align="right">(<a href="#readme-top">back to top</a>)</p> <!-- ROADMAP -->

## Roadmap

- [x] Automated background and shade corrections
- [x] Automated ROI detection and manual refinement
- [ ] Additional integration with alternative cell segmentation (Cellpose)
- [ ] Automated IRM-based deformation segmentation
- [ ] Add Changelog as updates roll in


See the [open issues](https://github.com/kscrudders/IT_SCAP/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p> <!-- CONTRIBUTING -->

Contributing: </p>
Contributions make this script more robust and easier to use. If you have suggestions:
* Fork the Project
* Create your Feature Branch (git checkout -b feature/YourFeature)
* Commit your Changes (git commit -m 'Added an awesome feature')
* Push to the Branch (git push origin feature/YourFeature)
* Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p> <!-- LICENSE -->

License: </p>
This project is distributed under GNU Genereal Public License. </p>
See LICENSE.txt for details.

<p align="right">(<a href="#readme-top">back to top</a>)</p> <!-- CONTACT -->
Contact </p>
Kevin Scrudders – kscrudders@gmail.com

Project Link: https://github.com/kscrudders/IT_SCAP

<p align="right">(<a href="#readme-top">back to top</a>)</p> <!-- ACKNOWLEDGMENTS -->

Acknowledgments
* The lab of Dr. Shalini T. Low-Nam
* The ever excellent MathWorks Documentation
* The code was developed from 2021-2025. During the later stages some code was drafted using ChatGPT. All code was reviewed, stress tested, and approved by me, Kevin.

<p align="right">(<a href="#readme-top">back to top</a>)</p>
