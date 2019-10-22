
A new method for classifying clast roundness based on the Radon transform

Moreno Chávez G.a, Jesús Villa b*, Sarocchi D.c, and Efrén González-Ramírezb

a Maestría en Ciencias del Procesamiento de la Información, Universidad Autónoma de Zacatecas, Av. Ramón López Velarde 801, Zacatecas 98000, México.
 b Unidad Académica de Ingeniería Eléctrica, Universidad Autónoma de Zacatecas, Av. Ramón López Velarde 801, Zacatecas 98000, México.
c Instituto de Geología / Fac. de Ingeniería UASLP, Av. Dr. M. Nava No 5, Zona Universitaria 78290, San Luis Potosí, México. 
*Corresponding author. Email address: jvillah@uaz.edu.mx  

_______________________________________________________________________
Abstract
In this paper, a new algorithm for clast’s roundness classification based on
 the Radon Transform is presented. The degree of roundness is determined by 
processing the sinogram of clasts images. The algorithm consists in applying 
two low pass filters to the sinogram, obtaining the inverse Radon transform and 
comparing the filtered image with the original one. For rounded particles, 
the difference between the original image and the filtered image will be small. 
On the other hand, considering angular clasts, the difference will be greater 
than in the case of the rounded ones, due to the presence of high-frequency components.
The comparison consists of the subtraction of the original image and the filtered image. 
Since the images are binary, the difference is an image with topologically unconnected 
regions that corresponds to the particle’s edges. The percentage difference between the 
original and the processed image and the number of regions are used to classify the 
clasts morphology. The results have been tested using comparison charts designed for 
visual roundness estimation (Powers, 1982). Two cutting frequencies, one to classify 
well-rounded, rounded and sub-rounded clasts and another one for angular and sub-angular
classes have been used. The proposed algorithm performs well in distinguish the roundness 
classes of the Powers’s visual charts. Similarly, the results provided by the algorithm have 
been compared with the classification performed by an expert. The algorithm attributes 92% of
the clasts to the same classes indicated by the expert. In this paper, we propose Gaussian 
models useful to classify the particles based on the Powers’s classes are also proposed.
User-Friendly software that allows clasts morphology classification by means of the Radon
algorithm, has been developed in the MATLAB platform. This software can be freely downloaded
from our web pages.


_______________________________________________________________________
We have developed a software, called RadonS, 
which contains the code of the developed method 
and the graphical interface. The code has been 
programmed in MATLAB 2014b. The main file is RadonShape.m. 
This file contains multiple functions. 
In the RadonShape_OpeningFcn function the global 
variables and elements of the graphic interface are declared. 
The load_Callback function loads and binarizes the image. 
The Radon_Callback function is the main function and performs 
the method described in this manuscript. The Radon_setupmodels_Callback
and Radon_plotmodels_Callback functions adjust the means and variances 
of the Gaussian models and are plotted, respectively. The class_Callback 
function classifies the particles in the analyzed images.  
The functions that perform the export of images and data are expdata_Callback
 and expima_Callback. The rest of the functions are related to the graphic 
elements of the graphic interface. In addition to the main RadonShape.m file, 
we have incorporated the RadonShape.fig file which is used to modify the interface 
from the MATLAB GUI, and models, lpfilter and normarea2 which are called by the main file. 
_______________________________________________________________________
To avoid incompatibility problems it is suggested 
to use the same version of matlab. However, 
if there is a compatibility error, 
you can use the executable version or contact us 
by the mail: gamalielmch@gmail.com 
The images used in the manuscript are attached 
The user manual is also attached
_______________________________________________________________________
