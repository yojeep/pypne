#include "inputData.h"
#include "blockNet.h"
// #include <format>

inline double randomG() /// to randomly distribute the shape factors, in case of errors
{
    double x1, x2, w, y;
    do
    {
        x1 = 2. * rand() / RAND_MAX - 1.;
        x2 = 2. * rand() / RAND_MAX - 1.;
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.);
    w = sqrt((-2. * log(w)) / w);
    y = 0.00625 * (x1 * w + 5.);
    if (y > 0.049)
        y = 0.0625;
    return y;
}

#pragma once
struct PoreData
{
    std::vector<int> pore_index;
    std::vector<float> pore_x;
    std::vector<float> pore_y;
    std::vector<float> pore_z;
    std::vector<int> pore_connection_number;
    std::vector<float> pore_volume; // 修正拼写错误 volumne -> volume
    std::vector<float> pore_radius;
    std::vector<float> pore_shape_factor;
};

#pragma once
struct ThroatData
{
    std::vector<int> throat_index;
    std::vector<int> throat_pore1_index;
    std::vector<int> throat_pore2_index;
    std::vector<float> throat_radius;
    std::vector<float> throat_shape_factor;
    std::vector<float> throat_total_length;
    std::vector<float> throat_pore1_length;
    std::vector<float> throat_pore2_length;
    std::vector<float> throat_length;
    std::vector<float> throat_volume; // 修正拼写
};

// 顶层数据结构
#pragma once
struct PYPNEData
{
    PoreData pore;
    ThroatData throat;
};

inline PYPNEData createPYPNEData(int n_nodes, int n_throats)
{
    PYPNEData data;
    data.pore.pore_index.resize(n_nodes);
    data.pore.pore_x.resize(n_nodes);
    data.pore.pore_y.resize(n_nodes);
    data.pore.pore_z.resize(n_nodes);
    data.pore.pore_connection_number.resize(n_nodes);
    data.pore.pore_volume.resize(n_nodes);
    data.pore.pore_radius.resize(n_nodes);
    data.pore.pore_shape_factor.resize(n_nodes);

    data.throat.throat_index.resize(n_throats);
    data.throat.throat_pore1_index.resize(n_throats);
    data.throat.throat_pore2_index.resize(n_throats);
    data.throat.throat_radius.resize(n_throats);
    data.throat.throat_shape_factor.resize(n_throats);
    data.throat.throat_total_length.resize(n_throats);
    data.throat.throat_pore1_length.resize(n_throats);
    data.throat.throat_pore2_length.resize(n_throats);
    data.throat.throat_length.resize(n_throats);
    data.throat.throat_volume.resize(n_throats);

    return data;
}
template <typename T>
class voxelImageT_PYPNE : public voxelImageT<T>
{
public:
    voxelImageT_PYPNE() = default;
    void setData(std::vector<T> &data)
    {
        voxelImageT<T>::data_ = std::move(data);
    }
    void setNnn(int3 nnn)
    {
        voxelImageT<T>::nnn_ = nnn;
    }
    void setNij(int nij)
    {
        voxelImageT<T>::nij_ = nij;
    }
};

class inputDataNE_PYPNE : public inputDataNE
{
public:
    inputDataNE_PYPNE() = default; // 或 inputData_PYPNE() {}

    // 带参数的构造函数
    inputDataNE_PYPNE(const int nx, const int ny, const int nz, double resolution, voxelImageT_PYPNE<unsigned char> &vm, std::string name)
        : inputDataNE()
    {
        this->nx = nx;
        this->ny = ny;
        this->nz = nz;
        this->vxlSize = resolution;
        this->VImage = std::move(vm);
        this->X0 = {0., 0., 0.};
        this->invalidSeg = {-10000, 255, 0};
        this->_rockTypes.push_back(poroRange("void", 0, 0));
        this->segValues.resize(256, this->_rockTypes.size());
        this->segValues[0] = 0;
        this->nBP6 = 2;
        this->name_ = name;
    }
    void setData(std::vector<std::pair<std::string, std::string>> data) { data_ = data; }
};

auto compute_PYPNEData_vector(blockNetwork &bn)
{

    /// pnflow uses the following indexes: [0:nB(=2)] for boundary nodes,
    /// [nB:nP+nB] for internal nodes and throat indices start afterwards.
    /// here, in Statoil format, all these are subtracted by 1 (starting from -1).

    /// ### First we compute the classical network model parameters,
    /// #### throat radii, shape factors and lengths,
    int nTrots = bn.nTrots;
    int nNodes = bn.nNodes;
    auto &throatIs = bn.throatIs;
    auto &poreIs = bn.poreIs;
    auto &cg = bn.cg;
    /// ### First we compute the classical network model parameters,
    /// #### throat radii, shape factors and lengths,
    vector<double> t_radiuss(nTrots, 0.);    //
    vector<double> t_shapeFacts(nTrots, 0.); //
    vector<double> t_lengthP1toP2s(nTrots, 0.);
    vector<double> t_lp1s(nTrots, 0.);  //
    vector<double> t_lp2s(nTrots, 0.);  //
    vector<double> t_ltrot(nTrots, 0.); // throat portion of t_lengthP1toP2s

    /// #### pore radii, shape factors and lengths.
    vector<double> p_radiuss(nNodes, 0.);
    vector<double> p_shape1s(nNodes, 0.);
    vector<double> p_physlengths(nNodes, 0.);

    dbl3 checkSumAt(0., 0., 0.);
    cout << "\ncalcThroats:" << endl;
    int lengthP1toP2Warnings = 0;
    double nBelowAllowedG(0.), nAboveAllowedG(0.), totalArea(0.);

    /// ### Compute throat  parameters
    for (int ti = 0; ti < nTrots; ++ti)
    {
        throatNE &tr = *throatIs[ti];

        double lthroat = 0;
        double lp1 = 0;
        double lp2 = 0;

        if (tr.surfaceArea == 0)
            tr.surfaceArea = 6;
        if (mag(tr.CrosArea) < 0.01)
            tr.CrosArea[0] = 0.1;
        checkSumAt += tr.CrosArea;

        /// - compute distance between the throat centre and each of the two adjacent pore centres
        double lpt1 = ((tr.e1 < 2) ? (tr.e1 == 0 ? tr.mb22()->fi : cg.nx - tr.mb22()->fi) : dist(poreIs[tr.e1]->mb, tr.mb22()));
        double rp1 = std::max((tr.e1 < 2) ? tr.mb22()->R : poreIs[tr.e1]->mb->R, 1.0f);
        double lpt2 = ((tr.e2 < 2) ? (tr.e2 == 0 ? tr.mb22()->fi : cg.nx - tr.mb22()->fi) : dist(poreIs[tr.e2]->mb, tr.mb22()));
        double rp2 = std::max((tr.e2 < 2) ? tr.mb22()->R : poreIs[tr.e2]->mb->R, 1.0f);

        /// - throat radius is the radius of the largest maximal sphere  on the throat surface
        double rr = std::max(tr.mb22()->R, 0.5f);
        rr = std::min(std::min(rr, rp1), rp2);
        t_radiuss[ti] = rr + 0.5 * (0.5 - double(rand()) / RAND_MAX);
        /// - throat total length is the sum of the two half-throat lengths
        double lengthP1toP2 = lpt1 + lpt2;
        if (lengthP1toP2 < 3.)
            lengthP1toP2 = 3.01, ++lengthP1toP2Warnings;

        /// - each pore is given 67% of the total throat length, the rest is called the throat elngth
        lp1 = lpt1 * 0.67;
        lp2 = lpt2 * 0.67;
        if (tr.e1 < 2)
            lp1 = 1;
        if (tr.e2 < 2)
            lp2 = 1;
        lthroat = lengthP1toP2 - lp1 - lp2;

        if (lthroat < 0.0000001)
            lthroat = 1;

        t_shapeFacts[ti] = rr * rr / 4. / mag(tr.CrosArea); ///- new throat shape factor definition G = R^2/4A

        if (t_shapeFacts[ti] >= 0.09)
        {
            t_shapeFacts[ti] = std::min(0.079, t_shapeFacts[ti] / 2.);
            nAboveAllowedG += mag(tr.CrosArea);
        } //. shape factor can not be this large, error: probably shared throat, temporary fix to handle in the flow code

        if (t_shapeFacts[ti] < 0.01) //. shape factor can not be this small,
        {
            t_shapeFacts[ti] = std::max(randomG(), 0.01);
            nBelowAllowedG += mag(tr.CrosArea);
        }

        totalArea += mag(tr.CrosArea);

        t_lengthP1toP2s[ti] = lengthP1toP2 * 1.;
        t_lp1s[ti] = lp1 * 1.;
        t_lp2s[ti] = lp2 * 1.;
        t_ltrot[ti] = lthroat * 1;
    }
    cout << " P1-to-P2 length < 3    for " << lengthP1toP2Warnings << " throats" << endl;
    cout << " shapefactor: belowAllowedG " << nBelowAllowedG / totalArea * 100 << "%   aboveAllowedG " << nAboveAllowedG / totalArea * 100 << "%" << endl;
    cout << " checkSumAt: " << checkSumAt << endl;

    cout << "calc Pores" << endl;

    /// ### Compute pore  parameters
    for (int pid = 2; pid < nNodes; ++pid)
    {
        poreNE &por = *poreIs[pid];
        double radius = por.mb->R;
        p_radiuss[pid] = radius;
        if (por.surfaceArea < 1)
            por.surfaceArea = 6;
        if (por.volumn < 1)
            por.volumn = 1;

        /// - pore shape factor is computed from a weighted average of its throat shape factors
        double shapeFactor(5.e-38), SumTArea(1e-36);
        for (const auto &bi : por.contacts)
        {
            throatNE &tr = *throatIs[bi.second];
            shapeFactor += t_shapeFacts[bi.second] * mag(tr.CrosArea);
            SumTArea += mag(tr.CrosArea);
        }
        shapeFactor /= SumTArea;

        double porA(radius * radius / 4. / shapeFactor);
        SumTArea += porA;
        double pVol = por.volumn;
        por.volumn = pVol * porA / SumTArea;

        for (const auto &bi : por.contacts)
        {
            throatNE &tr = *throatIs[bi.second];
            tr.volumn += pVol * mag(tr.CrosArea) / SumTArea;
        }

        p_shape1s[pid] = shapeFactor; //(por.volumn*poreLengthMax*2/por.surfaceArea/por.surfaceArea);
    }

    const double dx = cg.vxlSize;
    int n_throats = int(throatIs.size());
    int n_nodes = int(poreIs.size() - 2);
    PYPNEData pypnedata = createPYPNEData(n_nodes, n_throats);
    for (int ti = 0; ti < nTrots; ++ti)
    {
        const throatNE &tr = *throatIs[ti];
        pypnedata.throat.throat_index[ti] = ti + 1;
        pypnedata.throat.throat_pore1_index[ti] = int(tr.e1 - 1);
        pypnedata.throat.throat_pore2_index[ti] = int(tr.e2 - 1);
        pypnedata.throat.throat_radius[ti] = tr.radius() * dx;
        pypnedata.throat.throat_shape_factor[ti] = t_shapeFacts[ti];
        pypnedata.throat.throat_total_length[ti] = t_lengthP1toP2s[ti] * dx;
        pypnedata.throat.throat_pore1_length[ti] = t_lp1s[ti] * dx;
        pypnedata.throat.throat_pore2_length[ti] = t_lp2s[ti] * dx;
        pypnedata.throat.throat_length[ti] = t_ltrot[ti] * dx;
        pypnedata.throat.throat_volume[ti] = tr.volumn * dx * dx * dx;
    }

    for (int pid = 2; pid < nNodes; ++pid) //. 0th and first are inlet and outlet elements
    {
        const poreNE &por = *poreIs[pid];
        pypnedata.pore.pore_index[pid - 2] = pid - 1;
        pypnedata.pore.pore_x[pid - 2] = por.mb->fi * dx;
        pypnedata.pore.pore_y[pid - 2] = por.mb->fj * dx;
        pypnedata.pore.pore_z[pid - 2] = por.mb->fk * dx;
        pypnedata.pore.pore_connection_number[pid - 2] = int(por.contacts.size());
        pypnedata.pore.pore_volume[pid - 2] = por.volumn * dx * dx * dx;
        pypnedata.pore.pore_radius[pid - 2] = por.radius() * dx;
        pypnedata.pore.pore_shape_factor[pid - 2] = p_shape1s[pid];
    }

    // /// care entries
    // std::vector<int> result;
    // int inlet = 0, outlet = 0;
    // for (std::map<int, int>::const_iterator bi = por.contacts.begin(); bi != por.contacts.end(); ++bi) {
    //     const throatNE* tb = throatIs[bi->second];
    //     if (tb->e1 == pid) {
    //         if (tb->e2 == 0)
    //             inlet = 1;
    //         else if (tb->e2 == 1)
    //             outlet = 1;
    //         result.push_back(tb->e2 - 1);
    //     } else {
    //         if (tb->e1 == 0)
    //             inlet = 1;
    //         else if (tb->e1 == 1)
    //             outlet = 1;
    //         result.push_back(tb->e1 - 1);
    //     }
    // }

    // // 第二部分：存储入口/出口标志
    // result.push_back(inlet);
    // result.push_back(outlet);

    // // 第三部分：存储喉道原始索引
    // for (std::map<int, int>::const_iterator bi = por.contacts.begin(); bi != por.contacts.end(); ++bi) {
    //     result.push_back(bi->second + 1);
    // }
    return pypnedata;
}

#include <sstream>
#include <iomanip> // 用于格式控制（如 setw, scientific）

auto compute_PYPNEData_string(blockNetwork &bn)
{
    /// pnflow uses the following indexes: [0:nB(=2)] for boundary nodes,
    /// [nB:nP+nB] for internal nodes and throat indices start afterwards.
    /// here, in Statoil format, all these are subtracted by 1 (starting from -1).

    /// ### First we compute the classical network model parameters,
    /// #### throat radii, shape factors and lengths,
    int nTrots = bn.nTrots;
    int nNodes = bn.nNodes;
    auto &throatIs = bn.throatIs;
    auto &poreIs = bn.poreIs;
    auto &cg = bn.cg;
    /// ### First we compute the classical network model parameters,
    /// #### throat radii, shape factors and lengths,
    vector<double> t_radiuss(nTrots, 0.);    //
    vector<double> t_shapeFacts(nTrots, 0.); //
    vector<double> t_lengthP1toP2s(nTrots, 0.);
    vector<double> t_lp1s(nTrots, 0.);  //
    vector<double> t_lp2s(nTrots, 0.);  //
    vector<double> t_ltrot(nTrots, 0.); // throat portion of t_lengthP1toP2s

    /// #### pore radii, shape factors and lengths.
    vector<double> p_radiuss(nNodes, 0.);
    vector<double> p_shape1s(nNodes, 0.);
    vector<double> p_physlengths(nNodes, 0.);

    dbl3 checkSumAt(0., 0., 0.);
    cout << "\ncalcThroats:" << endl;
    int lengthP1toP2Warnings = 0;
    double nBelowAllowedG(0.), nAboveAllowedG(0.), totalArea(0.);

    /// ### Compute throat  parameters
    for (int ti = 0; ti < nTrots; ++ti)
    {
        throatNE &tr = *throatIs[ti];

        double lthroat = 0;
        double lp1 = 0;
        double lp2 = 0;

        if (tr.surfaceArea == 0)
            tr.surfaceArea = 6;
        if (mag(tr.CrosArea) < 0.01)
            tr.CrosArea[0] = 0.1;
        checkSumAt += tr.CrosArea;

        /// - compute distance between the throat centre and each of the two adjacent pore centres
        double lpt1 = ((tr.e1 < 2) ? (tr.e1 == 0 ? tr.mb22()->fi : cg.nx - tr.mb22()->fi) : dist(poreIs[tr.e1]->mb, tr.mb22()));
        double rp1 = std::max((tr.e1 < 2) ? tr.mb22()->R : poreIs[tr.e1]->mb->R, 1.0f);
        double lpt2 = ((tr.e2 < 2) ? (tr.e2 == 0 ? tr.mb22()->fi : cg.nx - tr.mb22()->fi) : dist(poreIs[tr.e2]->mb, tr.mb22()));
        double rp2 = std::max((tr.e2 < 2) ? tr.mb22()->R : poreIs[tr.e2]->mb->R, 1.0f);

        /// - throat radius is the radius of the largest maximal sphere  on the throat surface
        double rr = std::max(tr.mb22()->R, 0.5f);
        rr = std::min(std::min(rr, rp1), rp2);
        t_radiuss[ti] = rr + 0.5 * (0.5 - double(rand()) / RAND_MAX);
        /// - throat total length is the sum of the two half-throat lengths
        double lengthP1toP2 = lpt1 + lpt2;
        if (lengthP1toP2 < 3.)
            lengthP1toP2 = 3.01, ++lengthP1toP2Warnings;

        /// - each pore is given 67% of the total throat length, the rest is called the throat elngth
        lp1 = lpt1 * 0.67;
        lp2 = lpt2 * 0.67;
        if (tr.e1 < 2)
            lp1 = 1;
        if (tr.e2 < 2)
            lp2 = 1;
        lthroat = lengthP1toP2 - lp1 - lp2;

        if (lthroat < 0.0000001)
            lthroat = 1;

        t_shapeFacts[ti] = rr * rr / 4. / mag(tr.CrosArea); ///- new throat shape factor definition G = R^2/4A

        if (t_shapeFacts[ti] >= 0.09)
        {
            t_shapeFacts[ti] = std::min(0.079, t_shapeFacts[ti] / 2.);
            nAboveAllowedG += mag(tr.CrosArea);
        } //. shape factor can not be this large, error: probably shared throat, temporary fix to handle in the flow code

        if (t_shapeFacts[ti] < 0.01) //. shape factor can not be this small,
        {
            t_shapeFacts[ti] = std::max(randomG(), 0.01);
            nBelowAllowedG += mag(tr.CrosArea);
        }

        totalArea += mag(tr.CrosArea);

        t_lengthP1toP2s[ti] = lengthP1toP2 * 1.;
        t_lp1s[ti] = lp1 * 1.;
        t_lp2s[ti] = lp2 * 1.;
        t_ltrot[ti] = lthroat * 1;
    }
    cout << " P1-to-P2 length < 3    for " << lengthP1toP2Warnings << " throats" << endl;
    cout << " shapefactor: belowAllowedG " << nBelowAllowedG / totalArea * 100 << "%   aboveAllowedG " << nAboveAllowedG / totalArea * 100 << "%" << endl;
    cout << " checkSumAt: " << checkSumAt << endl;

    cout << "calc Pores" << endl;

    /// ### Compute pore  parameters
    for (int pid = 2; pid < nNodes; ++pid)
    {
        poreNE &por = *poreIs[pid];
        double radius = por.mb->R;
        p_radiuss[pid] = radius;
        if (por.surfaceArea < 1)
            por.surfaceArea = 6;
        if (por.volumn < 1)
            por.volumn = 1;

        /// - pore shape factor is computed from a weighted average of its throat shape factors
        double shapeFactor(5.e-38), SumTArea(1e-36);
        for (const auto &bi : por.contacts)
        {
            throatNE &tr = *throatIs[bi.second];
            shapeFactor += t_shapeFacts[bi.second] * mag(tr.CrosArea);
            SumTArea += mag(tr.CrosArea);
        }
        shapeFactor /= SumTArea;

        double porA(radius * radius / 4. / shapeFactor);
        SumTArea += porA;
        double pVol = por.volumn;
        por.volumn = pVol * porA / SumTArea;

        for (const auto &bi : por.contacts)
        {
            throatNE &tr = *throatIs[bi.second];
            tr.volumn += pVol * mag(tr.CrosArea) / SumTArea;
        }

        p_shape1s[pid] = shapeFactor; //(por.volumn*poreLengthMax*2/por.surfaceArea/por.surfaceArea);
    }

    const double dx = cg.vxlSize;
    cout << "Writing throats";
    cout.flush();

    // 保存所有生成内容的容器
    std::string link1Content, link2Content, node1Content, node2Content;

    { /// ### 生成 _link1.dat 的内容到字符串
        std::stringstream ss;
        ss << std::setw(6) << throatIs.size() << "\n";
        for (int ti = 0; ti < int(throatIs.size()); ++ti)
        {
            const throatNE &tr = *throatIs[ti];
            ss << std::setw(6) << ti + 1 << " "
               << std::setw(6) << int(tr.e1 - 1) << " "
               << std::setw(6) << int(tr.e2 - 1) << " "
               << std::scientific << tr.radius() * dx << " "
               << std::scientific << t_shapeFacts[ti] << " "
               << std::scientific << t_lengthP1toP2s[ti] * dx << "\n";
        }
        link1Content = ss.str();
    }
    // { /// ### 生成 _link1.dat 的内容到字符串
    //     link1Content = std::format("{:6d}\n", throatIs.size());
    //     for (int ti = 0; ti < int(throatIs.size()); ++ti)
    //     {
    //         const throatNE &tr = *throatIs[ti];
    //         link1Content += std::format(
    //             "{:6d} {:6d} {:6d} {:.6e} {:.6e} {:.6e}\n",
    //             ti + 1,
    //             int(tr.e1 - 1),
    //             int(tr.e2 - 1),
    //             tr.radius() * dx,
    //             t_shapeFacts[ti],
    //             t_lengthP1toP2s[ti] * dx);
    //     }
    // }

    { /// ### 生成 _link2.dat 的内容到字符串
        std::stringstream ss;
        for (int ti = 0; ti < int(throatIs.size()); ++ti)
        {
            const throatNE &tr = *throatIs[ti];
            ss << std::setw(6) << ti + 1 << " "
               << std::setw(6) << int(tr.e1 - 1) << " "
               << std::setw(6) << int(tr.e2 - 1) << " "
               << std::scientific << t_lp1s[ti] * dx << " "
               << std::scientific << t_lp2s[ti] * dx << " "
               << std::scientific << t_ltrot[ti] * dx << " "
               << std::scientific << tr.volumn * dx * dx * dx << " "
               << 0. << "\n";
        }
        link2Content = ss.str();
    }

    // { /// ### 生成 _link2.dat 的内容到字符串
    //     for (int ti = 0; ti < int(throatIs.size()); ++ti)
    //     {
    //         const throatNE &tr = *throatIs[ti];
    //         link2Content += std::format(
    //             "{:6d} {:6d} {:6d} {:.6e} {:.6e} {:.6e} {:.6e} {:.1f}\n",
    //             ti + 1,
    //             int(tr.e1 - 1),
    //             int(tr.e2 - 1),
    //             t_lp1s[ti] * dx,
    //             t_lp2s[ti] * dx,
    //             t_ltrot[ti] * dx,
    //             tr.volumn * dx * dx * dx,
    //             0.0);
    //     }
    // }

    cout << ".\n";
    cout.flush();

    cout << "Writing pores";
    cout.flush();

    { /// ### write  _node1.dat file 到字符串
        std::stringstream ss;
        // 写入头部（总节点数、nx, ny, nz）
        ss << std::setw(6) << (poreIs.size() - 2) << "  "
           << std::scientific << cg.nx * dx << "  "
           << std::scientific << cg.ny * dx << "  "
           << std::scientific << cg.nz * dx << "\n";

        for (int pid = 2; pid < int(poreIs.size()); ++pid) // 0和1是入口和出口
        {
            const poreNE &por = *poreIs[pid];

            // 写入节点ID、坐标、contact数量
            ss << std::setw(6) << int(pid - 1) << " "
               << std::scientific << (por.mb->fi * dx + 0. * cg.X0[0]) << " "
               << std::scientific << (por.mb->fj * dx + 0. * cg.X0[1]) << " "
               << std::scientific << (por.mb->fk * dx + 0. * cg.X0[2]) << " "
               << std::setw(3) << int(por.contacts.size());

            int inlet = 0, outlet = 0;
            // 遍历每个contact，写入连接的节点ID
            for (auto bi = por.contacts.begin(); bi != por.contacts.end(); ++bi)
            {
                const throatNE *tb = throatIs[bi->second];
                if (tb->e1 == pid)
                {
                    // 检查是否连接到入口（0）或出口（1）
                    if (tb->e2 == 0)
                        inlet = 1;
                    else if (tb->e2 == 1)
                        outlet = 1;

                    // 写入连接的节点ID（tb->e2 -1）
                    ss << "\t" << std::setw(5) << int(tb->e2 - 1);
                }
                else
                {
                    // 检查反向连接
                    if (tb->e1 == 0)
                        inlet = 1;
                    else if (tb->e1 == 1)
                        outlet = 1;

                    // 写入连接的节点ID（tb->e1 -1）
                    ss << "\t" << std::setw(5) << int(tb->e1 - 1);
                }
            }

            // 写入入口/出口标记（占位符）
            ss << "\t" << std::setw(5) << inlet
               << "\t" << std::setw(5) << outlet;

            // 再次遍历contacts，写入throat的ID（bi->second +1）
            for (auto bi = por.contacts.begin(); bi != por.contacts.end(); ++bi)
            {
                ss << "\t" << std::setw(5) << int(bi->second + 1);
            }

            // 换行
            ss << "\n";
        }
        node1Content = ss.str();
    }
    // { /// ### write _node1.dat file 到字符串
    //     // 写入头部（总节点数、nx, ny, nz）
    //     node1Content += std::format("{:6d}  {:.6e}  {:.6e}  {:.6e}\n",
    //                                 poreIs.size() - 2,
    //                                 cg.nx * dx,
    //                                 cg.ny * dx,
    //                                 cg.nz * dx);

    //     for (int pid = 2; pid < int(poreIs.size()); ++pid) // 0和1是入口和出口
    //     {
    //         const poreNE &por = *poreIs[pid];

    //         // 写入节点基本信息
    //         node1Content += std::format("{:6d} {:.6e} {:.6e} {:.6e} {:3d}",
    //                                     pid - 1,
    //                                     por.mb->fi * dx + 0. * cg.X0[0],
    //                                     por.mb->fj * dx + 0. * cg.X0[1],
    //                                     por.mb->fk * dx + 0. * cg.X0[2],
    //                                     int(por.contacts.size()));

    //         int inlet = 0, outlet = 0;
    //         // 第一次遍历contacts：写入连接的节点ID
    //         for (const auto &[_, throat_idx] : por.contacts)
    //         {
    //             const throatNE *tb = throatIs[throat_idx];
    //             const bool is_e1 = (tb->e1 == pid);
    //             const int other_node = is_e1 ? tb->e2 : tb->e1;

    //             // 检查入口/出口
    //             if (other_node == 0)
    //                 inlet = 1;
    //             else if (other_node == 1)
    //                 outlet = 1;

    //             // 写入连接节点ID
    //             node1Content += std::format("\t{:5d}", other_node - 1);
    //         }

    //         // 写入入口/出口标记
    //         node1Content += std::format("\t{:5d}\t{:5d}", inlet, outlet);

    //         // 第二次遍历contacts：写入throat ID
    //         for (const auto &[_, throat_idx] : por.contacts)
    //         {
    //             node1Content += std::format("\t{:5d}", throat_idx + 1);
    //         }

    //         node1Content += "\n";
    //     }
    // }

    { /// ### 生成 _node2.dat 的内容到字符串
        std::stringstream ss;
        for (int pid = 2; pid < int(poreIs.size()); ++pid)
        {
            const poreNE &por = *poreIs[pid];
            ss << std::setw(6) << int(pid - 1) << " "
               << std::scientific << por.volumn * dx * dx * dx << " "
               << std::scientific << por.radius() * dx << " "
               << std::scientific << p_shape1s[pid] << " "
               << 0. << "\n";
        }
        node2Content = ss.str();
    }

    // { /// ### 生成 _node2.dat 的内容到字符串
    //     for (int pid = 2; pid < int(poreIs.size()); ++pid)
    //     {
    //         const poreNE &por = *poreIs[pid];
    //         node2Content += std::format(
    //             "{:6d} {:.6e} {:.6e} {:.6e} {:.1f}\n",
    //             pid - 1,
    //             por.volumn * dx * dx * dx,
    //             por.radius() * dx,
    //             p_shape1s[pid],
    //             0.0);
    //     }
    // }

    cout << ".\n";
    cout.flush();
    return std::make_tuple(link1Content, link2Content, node1Content, node2Content);
}