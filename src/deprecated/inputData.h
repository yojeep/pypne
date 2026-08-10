// #ifndef CONFIG_H
// #define CONFIG_H

// #include "InputFile.h"
// #include "voxelImage.h"
// #include "ElementGNE.h"



// #ifndef PORORANGE_H
// #define PORORANGE_H
// class poroRange : public std::pair<unsigned char, unsigned char>
// {

// public:
// 	poroRange(unsigned char lower, unsigned char upper) : std::pair<unsigned char, unsigned char>(lower, upper) {};
// 	poroRange(std::string nam, unsigned char lower, unsigned char upper)
// 		: std::pair<unsigned char, unsigned char>(lower, upper), name(nam) {};
// 	poroRange(const poroRange &copyt) : std::pair<unsigned char, unsigned char>(copyt), name(copyt.name) {};
// 	bool outside(unsigned char value) { return this->first > value || value > this->second; };
// 	std::string name;
// };
// #endif // PORORANGE_H

// inline int createSample_input_nextract(std::string fnam, std::string opts)
// {
// 	if (fnam.empty())
// 		fnam = "vxlImage.mhd";

// 	ensure(!std::ifstream(fnam), "\n\nFile " + fnam + " exists,\n" + "  to run simulation: rerun with " + fnam + " as the only argument \n" + "  to regenerate: delete it and try again,\n or provide a different file name:\n" + "   pnextract -g input_pnextract.mhd\n", -1);

// 	std::ofstream of(fnam);
// 	if (opts == "-g")
// 	{
// 		of << "ObjectType =  Image\n"
// 		   << "NDims =       3\n"
// 		   << "ElementType = MET_UCHAR\n"
// 		   << "ElementByteOrderMSB = False\n"
// 		   << "ElementNumberOfChannels = 1\n"
// 		   << "CompressedData = True\n\n"
// 		   << "HeaderSize = 0\n"
// 		   << "DimSize =    	1000	1000	1000\n"
// 		   << "ElementSize = 	1.6 	1.6 	1.6\n"
// 		   << "Offset =      	0   	0   	0\n"
// 		   << "\n"
// 		   << "ElementDataFile = input_image.raw.gz\n"
// 		   << "\n\n"
// 		   << "//! The above keywords are compatible with mhd format and \n"
// 		   << "//! can be used to open the Image in Fiji/ImageJ or Paraview\n"
// 		   << "//! The following commands are optional, remove the \"//\" to activate them\n"
// 		   << "\n"
// 		   << "//DefaultImageFormat = tif\n"
// 		   << "\n"
// 		   << "//!______________  image processing  commands _________________\n"
// 		   << "\n"
// 		   << "//! crop image to  [ Nxyz_begin  Nxyz_end )"
// 		   << "//cropD                0 0 0    300 300 300 \n"
// 		   << "\n"
// 		   << "//! flip x direction with y or z \n"
// 		   << "//direction z\n"
// 		   << "\n"
// 		   << "//! manipulate voxel values:\n"
// 		   << "//!    range   [start...end] -> value \n"
// 		   << "//replaceRange   0     127       0 \n"
// 		   << "//replaceRange   128   255       1 \n"
// 		   << "\n"
// 		   << "//! threshold image: range -> 0 (void-space), rest->1 (solid) \n"
// 		   << "//threshold   0  128 \n"
// 		   << "\n\n\n";
// 		of << "//!_______________  network extraction keywords __________________\n"
// 		   << "//!______(should be after image processing  commands above) ______\n"
// 		   << "\n"
// 		   << "//title:   output_network"
// 		   << "\n"
// 		   << "//write_all:	true; // use `write_all` is a memorable alternative to all visualization keywords "
// 		   << "// write_radius:	true\n"
// 		   << "// write_statistics:	true\n"
// 		   << "// write_elements:	true\n"
// 		   << "// write_poreMaxBalls:	true\n"
// 		   << "// write_throatMaxBalls:	true\n"
// 		   << "// write_throats:	true\n"
// 		   //<<"// write_poroats:	true\n" // results in seg fault
// 		   << "// write_hierarchy:	true\n"
// 		   << "// write_medialSurface:	true\n"
// 		   << "// write_throatHierarchy:	true\n"
// 		   << "// write_vtkNetwork:	true\n"
// 		   << std::endl;
// 	}
// 	else
// 		std::cout << "Error unknown option (first argument)" << std::endl;

// 	of.close();

// 	std::cout << " file " << fnam << " generated, edit: set Image size, name etc, and rerun\n";
// 	return 0;
// }

// class inputDataNE : public InputFile
// {
// public:
// 	inputDataNE() = default;
// 	inputDataNE(const std::string &fnam)
// 		: InputFile(fnam, false), nx(0), ny(0), nz(0), vxlSize(1), X0(0., 0., 0.), invalidSeg{-10000, 255, 0}
// 	{
// 	}

// 	void init(bool verbos = true)
// 	{


// 		nBP6 = getOr("multiDir", false) ? 6 : 2;

// 		std::cout << " voxel indices:" << std::endl;
// 		_rockTypes.push_back(poroRange("void", 0, 0));
// 		std::cout << "  " << 0 << ": void voxels " << std::endl;

// 		std::istringstream iss;

// 		for_i(_rockTypes)
// 		{
// 			auto &rt = _rockTypes[i];
// 			if (giv(rt.name + "_range", iss))
// 			{
// 				int lower, upper;
// 				iss >> lower >> upper;
// 				rt.first = lower;
// 				rt.second = upper;
// 				if (rt.first > rt.second)
// 					std::cout << "  Wrong entries for keyword \"" << rt.name + "_range" << "\":  lower value is higher than upper value" << std::endl;
// 			}

// 			for (size_t j = rt.first; j <= rt.second; ++j)
// 				segValues[j] = i;

// 			std::cout << "  " << rt.name << " voxel values: [" << int(rt.first) << " " << int(rt.second) << "]" << std::endl;
// 		}

// 		std::cout << "  Voxel value indices:";
// 		for (size_t i = 0; i <= 12; ++i)
// 			std::cout << " " << segValues[i];
// 		std::cout << " ... " << std::endl;
// 	}

// 	void readImage()
// 	{

// 		std::string fnam(fileName());
// 		if (fnam.empty())
// 		{
// 			giv("ElementDataFile", fnam) || giv("read", fnam);
// 		}
// 		// Assert(giv("ElementDataFile") || kwrd("read").size()

// 		// std::string dataType("mhd");  giv("ElementType", dataType) || giv("fileType", dataType);
// 		std::cout << " Image file: " << fnam << std::endl;

// 		// VImage.reset(nx,ny,nz,255);
// 		// if (dataType == "microct")	VImage.readMicroCT(fnam);
// 		// else if (dataType == "binary")  {if (!VImage.readBin(fnam))   alert("could not read binary image!",-1); }
// 		// else if (dataType == "ascii")   {if (!VImage.readAscii(fnam)) alert("could not read ascii image!",-1); }
// 		// else  { VImage.reset(0,0,0,255);
// 		readConvertFromHeader(VImage, fnam); // dataType="mhd"; }

// 		const std::string &vxlkys = kwrd("VxlPro");
// 		if (vxlkys.size())
// 			vxlProcess(vxlkys, VImage, "GNE:VxlPro");

// 		// if(dataType=="mhd")
// 		//{
// 		vxlSize = VImage.dx().x;
// 		X0 = VImage.X0();
// 		//}
// 		nx = VImage.nx();
// 		ny = VImage.ny();
// 		nz = VImage.nz();
// 		ensure(nz, "no image read", 2);
// 		VImage.printInfo();

// 		nInside = (long long)(nx)*ny * nz;
// 		std::cout << " siz:" << VImage.size3() << " vxlSize:" << vxlSize << " X0:" << X0 << std::endl;

// 		int2 outrange;
// 		if (giv("outside_range", outrange))
// 		{
// 			forAllcp(VImage) if (outrange.a <= (*cp) && (*cp) <= outrange.b)-- nInside;
// 			std::cout << " outside_range: " << outrange << "; inside_fraction: " << nInside / (double(nx) * ny * nz) << std::endl;
// 		}
// 	}



// public:
// 	int nx, ny, nz;
// 	int nBP6;
// 	double vxlSize;
// 	dbl3 X0;
// 	std::string imgfrmt;
// 	long long nInside;

// 	std::string flowBaseDir;

// 	stvec<poroRange> _rockTypes;
// 	voxelImage VImage;

// };

// #endif
