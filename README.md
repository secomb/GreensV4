# Green’s function method version 4 (non-GPU version)

In the Green’s function method for simulating solute transport from a network of vessels to a finite volume of tissue, vessels and tissue are treated as distributions of sources of sinks. The method is described in: Secomb TW, Hsu R, Park EY, Dewhirst MW. 2004 Green's function methods for analysis of oxygen delivery to tissue by microvascular networks. Ann. Biomed. Eng 32, 1519-1529.

Notes

1. This method uses an “infinite-domain solution,” in which the network of vessels and the associated tissue domain are effectively embedded in an infinite space with zero net exchange of solute across the tissue domain boundaries. Such an approach has two advantages: it is applicable to tissue domains of arbitrary shape, and it avoids artifacts that can occur when specific boundary conditions, such as the no-flux condition, are imposed.

2. This method can be used for multiple reacting solutes. It allows for diffusible and non-diffusible solutes. Diffusible solutes may be specified as permeating vessel walls or non-permeating. Oxygen is handled differently from other solutes to account for binding by hemoglobin. Intravascular resistance to solute transport can be represented as a function of vessel diameter, as is the case for oxygen. For other solutes, the permeability of the vessel wall to the solute can be specified.

3. The network of vessels is approximated by a set of straight uniform segments located inside a reference cuboidal (box) shape. This cuboid must include all tissue regions associated with the given network. Each segment is then divided into sub-segments with equal length. The midpoint of each sub-segment represents a source. The tissue domain is divided into small subregions, each centered on a tissue node point. The tissue node points form a three dimension matrix. Inside each subregion, the solute reaction rates are assumed to be uniform, and depend on the solute levels at the tissue node.

4. The tissue node points within the tissue domain are determined by the function outboun.cpp, which is called from analyzenet.cpp with argument 1 or 2 to specify the method. Method 1 finds the smallest convex region inside the cuboid which excludes tissue node points that have a distance to the nearest vessel greater than a value specified by the user (Lb in network.dat). Method 2 finds all tissue points within a distance lb of the vessels. The resulting region is the simulated tissue domain, which is embedded in an infinite domain. The whole cuboidal region can be included in the tissue domain by making Lb large enough.

5. We supply several sample problems in the sub-directories. For input file formatting, please see the files in these directories. These files contain further explanatory comments.

6. Input files must be placed in the same working directory where the C code resides. They can be copied from the folders containing sample problems. The following input files are used:  
Network.dat This file specifies the network structure, the flow rate and hematocrit of all segments and the PO2 and solute concentrations of boundary nodes. (Note that these solute concentrations are used only for inflow nodes.)  
SoluteParams.dat This file gives solute transport and reaction parameters. Oxygen is handled differently than other solutes to account for binding by hemoglobin. Other solutes are assumed to be carried in solution and not to bind to blood components. This file contains oxygen transport parameters, which are needed only if the “oxygen” parameter is set to 1.  
ContourParams.dat At each time point, the program generates a postscript file showing the vessel network projected onto a single plane and the solute contours on that plane. If there are multiple solutes, the plots are generated for each plot as separate pages within the postscript file. The position and orientation of the plane are specified by ContourParams.dat, giving the coordinates of three corners of the plane. Also, the contour levels are specified in this file.  
IntravascRes.dat Specifies intravascular resistance to solute transport from blood to tissue. For oxygen, this depends on the vessel diameter. For other solutes, a wall permeability is specified by the file.  
VaryParams.dat [New in Version 4] The program provides the capability to specify multiple runs with different parameter values. In version 4, an optional VaryParams.dat file can be included. The following parameters can be varied:  
q0fac, diff[isp], tissparam[i][isp] - specify actual values  
solutefac[isp], intravascfac[isp] - specify multiplicative factors  
A run is performed for each specified set of parameters. An example is included in the “Brain” folder. Note that the values specified in VaryParams.dat override the values given in the other input data files.  
tissrate.cpp.dat This is a segment of C code supplied by the user that gives the rates of production mtiss[isp] of solute number isp in the tissue as a function of the solute concentrations c[1], c[2], …, c[nsp]. Also, the derivatives mptiss[isp] of mtiss[isp] with respect to c[isp] must be specified in this code. Caution: If this file is modified, the project must be rebuilt otherwise the old version will remain in the compiled code.  
postgreens.cpp.dat [New in Version 4] This user-supplied c code specifies computation of dependent variables from the results of greens. For example, concentrations of cytotoxic drugs can be used to compute cellular survival fractions at each point in the tissue.  
PostGreensParams.dat [New in Version 4] This file specifies parameters that are used in postgreens.cpp.  

7. The program generates several output files giving results in numerical form. It generates a folder called "Current" in the working folder (if it does not already exist), and generates postscript files showing the vessel network and the contours of each solute, in the "Current" folder. The contouring routine produces a color-shaded contour plot. It includes procedures to consolidate polygon of the same color, so as to reduce the size of the resulting output postscript files. The network is superimposed on the contour plot, and is plotted in order of increasing z with vessels outlined to give a pseudo three-dimensional effect. It can be viewed using a postscript viewer, such as GSview https://www.ghostscript.com/download/gsdnld.html.   Other output files: summary.out - Summary of runs; useful when VaryParams.dat is used.
VesselSources.out, VesselLevels.out, TissueSources.out, TissueLevels.out - Data on solute levels and source strengths at vessel points and at tissue points.
greens.exelem, greens.exenode, greens.com.txt [New in Version 4] - Files used to obtain 3D visualizations of the network using CMGUI. When CMGUI is started, use File - Open - com file, select greens.com.txt and hit “All” for the visualization. For details of cmgui, see: http://sourceforge.net/projects/cmiss/files/cmgui/cmgui-wx-2.8.0/

8. We have tested this package using Microsoft Visual C++ 2017 under Windows 7. For instructions on setting up a new project with this compiler, please see NewProject.txt. For error reporting and suggestions please contact Dr. Timothy W. Secomb, (520) 626-4513, email secomb@u.arizona.edu. We welcome your comments and suggestions.

9. This program is freely available for non-commercial use, provided appropriate acknowledgement is given. Commercial users please contact us before using this program. No assurance is given that it is free of errors and any use is at the user’s risk.

Version 4.0 changes

1. Method to vary blood flow rates throughout network. In versions 3 and below, the parameter q0 was included (SoluteParams.dat, line 10) to give the total inflow to the network. In version 4, a parameter q0fac is defined (SoluteParams.dat, line 10). When the program is run, all segment flows specified in Network.dat are multiplied by q0fac. This change is NOT compatible with SoluteParams.dat files used for versions 3 and below.
2. Number of tissue parameters per solute. In versions 3 and below, three tissue parameters at most could be defined per solute. In version 4, the number of parameters is defined as the second entry in line 11 of SoluteParams.dat. This change is NOT compatible with SoluteParams.dat files used for versions 3 and below.
3. Capability to specify multiple runs with different parameter values has been added. See above.
4. Capability for 3D visualization of networks using CMGUI has been added. See above.
5. Capability has been added to compute dependent quantities after greens has run. See above.
6. Summary output file from multiple runs (summary.out) has been added.
7. Output files are formatted differently. TissueLevels.dat, TissueSources.dat, VesselLevels.dat, VesselSources.dat now have one line per tissue point or vessel point, for ease of importing into spreadsheets.
8. Previous combined output file (GreensRes.out) is no longer generated.
9. The capability to restart a run has been removed (was not functional).
10. All output files and a copy of input files are placed in a folder called "Current."
11. The algorithm for non-diffusible solutes has been improved by using reaction rates based on updated levels of other solutes at each iteration of Newton method.
12. Error bounds have been modified to use relative errors instead of absolute errors.

Updated 15 June 2019