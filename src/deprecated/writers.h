// #pragma once
// /*---------------------------------------------------------------------------*\
// written by:
//         Ali Q Raeini  email: a.q.raeini@imperial.ac.uk  and
//         Tom Bultreys
// \*---------------------------------------------------------------------------*/

// // #include "blockNet.h"
// // #include "medialSurf.h"
// // #include "voxelImage.h"
// #include "ElementGNE.h"
// #include "globals.h"
// #include "typses.h"

// class medialSurface;
// class blockNetwork;

// TensorXXXDf32 ballRadiiToVoxel(const blockNetwork &mpn);
// TensorXXXDi32 VElemsPlusThroats(const blockNetwork &mpn);

// ///- written by Tom Bultreys:
// TensorXXXDi32 VThroats(const blockNetwork &mpn);
// TensorXXXDi32 VThroats(const blockNetwork &mpn, int firstSlice, int lastSlice);
// TensorXXXDi32 poreMaxBalls(const blockNetwork &mpn);
// TensorXXXDi32 poreMaxBalls(const blockNetwork &mpn, int firstSlice,
//                            int lastSlice);
// TensorXXXDi32 throatMaxBalls(const blockNetwork &mpn);
// TensorXXXDi32 throatMaxBalls(const blockNetwork &mpn, int firstSlice,
//                              int lastSlice);

// void vtuWriteMbMbs(std::string bNam, const std::vector<medialBall> &ballSpace,
//                    const std::vector<poreNE> &poreIs,
//                    const TensorXXXDi32 &VElems, double dx, dbl3 X0);
// void vtuWriteThroatMbMbs(std::string baseNam,
//                          const std::vector<throatNE> &throatIs,
//                          const std::vector<poreNE> &poreIs,
//                          const TensorXXXDi32 &VElems, double dx, dbl3 X0);

// void vtuWritePores(std::string bNam, const std::vector<poreNE> &poreIs,
//                    const std::vector<throatNE> &throatIs, double dx, dbl3 X0);
// void vtuWriteTHroatSpheres(std::string bNam, const std::vector<poreNE> &poreIs,
//                            const std::vector<throatNE> &throatIs, double dx,
//                            dbl3 X0);
// void vtuWriteThroats(std::string bNam, const std::vector<poreNE> &poreIs,
//                      const std::vector<throatNE> &throatIs, double dx, dbl3 X0);
