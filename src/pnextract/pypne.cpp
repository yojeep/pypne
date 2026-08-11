#include "config.hpp"
#include "pypne.hpp"
#include <string>
#include "types.hpp"

std::tuple<std::string, std::string, std::string, std::string>
get_network(blockNetwork &bn) {
  /// pnflow uses the following indexes: [0:nB(=2)] for boundary nodes,
  /// [nB:nP+nB] for internal nodes and throat indices start afterwards.
  /// here, in Statoil format, all these are subtracted by 1 (starting from -1).

  /// ### First we compute the classical network model parameters,
  /// #### throat radii, shape factors and lengths,
  int nTrots = bn.nTrots;
  int nNodes = bn.nNodes;
  auto &throatIs = bn.throatIs;
  auto &poreIs = bn.poreIs;
  auto &cg = bn.cfg;
  /// ### First we compute the classical network model parameters,
  /// #### throat radii, shape factors and lengths,
  std::vector<double> t_radiuss(nTrots, 0.);    //
  std::vector<double> t_shapeFacts(nTrots, 0.); //
  std::vector<double> t_lengthP1toP2s(nTrots, 0.);
  std::vector<double> t_lp1s(nTrots, 0.);  //
  std::vector<double> t_lp2s(nTrots, 0.);  //
  std::vector<double> t_ltrot(nTrots, 0.); // throat portion of t_lengthP1toP2s

  /// #### pore radii, shape factors and lengths.
  std::vector<double> p_radiuss(nNodes, 0.);
  std::vector<double> p_shape1s(nNodes, 0.);
  std::vector<double> p_physlengths(nNodes, 0.);

  Vector3f32 checkSumAt(0., 0., 0.);
  std::cout << "\ncalcThroats:" << std::endl;
  int lengthP1toP2Warnings = 0;
  double nBelowAllowedG(0.), nAboveAllowedG(0.), totalArea(0.);

  /// ### Compute throat  parameters
  for (int ti = 0; ti < nTrots; ++ti) {
    throatNE &tr = throatIs[ti];

    double lthroat = 0;
    double lp1 = 0;
    double lp2 = 0;

    if (tr.surfaceArea == 0)
      tr.surfaceArea = 6;

    if (tr.CrosArea.norm() < 0.01f)
      tr.CrosArea[0] = 0.1f;
    checkSumAt += tr.CrosArea;

    /// - compute distance between the throat centre and each of the two
    /// adjacent pore centres
    double lpt1 =
        ((tr.e1 < 2) ? (tr.e1 == 0 ? tr.mb22()->fx : cg.nx - tr.mb22()->fx)
                     : dist(poreIs[tr.e1].mb, tr.mb22()));
    double rp1 =
        std::max((tr.e1 < 2) ? tr.mb22()->R : poreIs[tr.e1].mb->R, 1.0f);
    double lpt2 =
        ((tr.e2 < 2) ? (tr.e2 == 0 ? tr.mb22()->fx : cg.nx - tr.mb22()->fx)
                     : dist(poreIs[tr.e2].mb, tr.mb22()));
    double rp2 =
        std::max((tr.e2 < 2) ? tr.mb22()->R : poreIs[tr.e2].mb->R, 1.0f);

    /// - throat radius is the radius of the largest maximal sphere  on the
    /// throat surface
    double rr = std::max(tr.mb22()->R, 0.5f);
    rr = std::min(std::min(rr, rp1), rp2);
    t_radiuss[ti] = rr + 0.5 * (0.5 - double(rand()) / RAND_MAX);
    /// - throat total length is the sum of the two half-throat lengths
    double lengthP1toP2 = lpt1 + lpt2;
    if (lengthP1toP2 < 3.)
      lengthP1toP2 = 3.01, ++lengthP1toP2Warnings;

    /// - each pore is given 67% of the total throat length, the rest is called
    /// the throat elngth
    lp1 = lpt1 * 0.67;
    lp2 = lpt2 * 0.67;
    if (tr.e1 < 2)
      lp1 = 1;
    if (tr.e2 < 2)
      lp2 = 1;
    lthroat = lengthP1toP2 - lp1 - lp2;

    if (lthroat < 0.0000001)
      lthroat = 1;

    t_shapeFacts[ti] =
        rr * rr / 4. /
        tr.CrosArea.norm(); ///- new throat shape factor definition G = R^2/4A

    if (t_shapeFacts[ti] >= 0.09) {
      t_shapeFacts[ti] = std::min(0.079, t_shapeFacts[ti] / 2.);
      nAboveAllowedG += tr.CrosArea.norm();
    } //. shape factor can not be this large, error: probably shared throat,
      // temporary fix to handle in the flow code

    if (t_shapeFacts[ti] < 0.01) //. shape factor can not be this small,
    {
      t_shapeFacts[ti] = randomG();
      nBelowAllowedG += tr.CrosArea.norm();
    }

    totalArea += tr.CrosArea.norm();

    t_lengthP1toP2s[ti] = lengthP1toP2 * 1.;
    t_lp1s[ti] = lp1 * 1.;
    t_lp2s[ti] = lp2 * 1.;
    t_ltrot[ti] = lthroat * 1;
  }
  std::cout << " P1-to-P2 length < 3    for " << lengthP1toP2Warnings
            << " throats" << std::endl;
  std::cout << " shapefactor: belowAllowedG "
            << nBelowAllowedG / totalArea * 100 << "%   aboveAllowedG "
            << nAboveAllowedG / totalArea * 100 << "%" << std::endl;
  std::cout << " checkSumAt: " << checkSumAt << std::endl;

  std::cout << "calc Pores" << std::endl;

  /// ### Compute pore  parameters
  for (int pid = 2; pid < nNodes; ++pid) {
    poreNE &por = poreIs[pid];
    double radius = por.mb->R;
    p_radiuss[pid] = radius;
    if (por.surfaceArea < 1)
      por.surfaceArea = 6;
    if (por.volumn < 1)
      por.volumn = 1;

    /// - pore shape factor is computed from a weighted average of its throat
    /// shape factors
    double shapeFactor(5.e-38), SumTArea(1e-36);
    for (const auto &bi : por.contacts) {
      throatNE &tr = throatIs[bi.second];
      shapeFactor += t_shapeFacts[bi.second] * tr.CrosArea.norm();
      SumTArea += tr.CrosArea.norm();
    }
    shapeFactor /= SumTArea;

    double porA(radius * radius / 4. / shapeFactor);
    SumTArea += porA;
    double pVol = por.volumn;
    por.volumn = pVol * porA / SumTArea;

    for (const auto &bi : por.contacts) {
      throatNE &tr = throatIs[bi.second];
      tr.volumn += pVol * tr.CrosArea.norm() / SumTArea;
    }

    p_shape1s[pid] =
        shapeFactor; //(por.volumn*poreLengthMax*2/por.surfaceArea/por.surfaceArea);
  }

  const double dx = cg.resolution;
  std::cout << "Writing throats";
  std::cout.flush();

  // 保存所有生成内容的容器
  std::string link1Content, link2Content, node1Content, node2Content;
  // int 6, float 15, space 1, return 1

  { // link1.dat
    std::format_to(std::back_inserter(link1Content), "{:6d}\n",
                   throatIs.size());

    for (int ti = 0; ti < int(throatIs.size()); ++ti) {
      const throatNE &tr = throatIs[ti];
      std::format_to(std::back_inserter(link1Content),
                     "{:6d} {:6d} {:6d} {:.6E} {:.6E} {:.6E}\n", ti + 1,
                     int(tr.e1 - 1), int(tr.e2 - 1), tr.radius() * dx,
                     t_shapeFacts[ti], t_lengthP1toP2s[ti] * dx);
    }
  }

  { // link2.dat
    for (int ti = 0; ti < int(throatIs.size()); ++ti) {
      const throatNE &tr = throatIs[ti];
      std::format_to(std::back_inserter(link2Content),
                     "{:6d} {:6d} {:6d} {:.6E} {:.6E} {:.6E} {:.6E} {:.6E}\n",
                     ti + 1, int(tr.e1 - 1), int(tr.e2 - 1), t_lp1s[ti] * dx,
                     t_lp2s[ti] * dx, t_ltrot[ti] * dx,
                     tr.volumn * dx * dx * dx, 0.0);
    }
  }

  std::cout << ".\n";
  std::cout.flush();

  std::cout << "Writing pores";
  std::cout.flush();
  { // node1.dat
    std::format_to(std::back_inserter(node1Content),
                   "{:6d} {:.6E} {:.6E} {:.6E}\n", poreIs.size() - 2,
                   cg.nx * dx, cg.ny * dx, cg.nz * dx);
    for (int pid = 2; pid < int(poreIs.size()); ++pid) {
      const poreNE &por = poreIs[pid];

      std::format_to(std::back_inserter(node1Content),
                     "{:6d} {:.6E} {:.6E} {:.6E} {:3d}", pid - 1,
                     por.mb->fx * dx, por.mb->fy * dx, por.mb->fz * dx,
                     int(por.contacts.size()));

      int inlet = 0, outlet = 0;
      for (auto bi = por.contacts.begin(); bi != por.contacts.end(); ++bi) {
        const throatNE &tb = throatIs[bi->second];

        int other = (tb.e1 == pid) ? tb.e2 : tb.e1;

        inlet |= (other == 0);
        outlet |= (other == 1);
        int neighbor = static_cast<int>(other - 1);

        std::format_to(std::back_inserter(node1Content), " {:5d}", neighbor);
      }

      std::format_to(std::back_inserter(node1Content), " {:5d} {:5d}", inlet,
                     outlet);

      for (auto bi = por.contacts.begin(); bi != por.contacts.end(); ++bi) {
        std::format_to(std::back_inserter(node1Content), " {:5d}",
                       int(bi->second + 1));
      }
      node1Content += '\n';
    }
  }

  { // node2.dat
    for (int pid = 2; pid < int(poreIs.size()); ++pid) {
      const poreNE &por = poreIs[pid];
      std::format_to(std::back_inserter(node2Content),
                     "{:6d} {:.6E} {:.6E} {:.6E} {:.6E}\n", pid - 1,
                     por.volumn * dx * dx * dx, por.radius() * dx,
                     p_shape1s[pid], 0.);
    }
  }

  std::cout << ".\n";
  std::cout.flush();
  return std::make_tuple(link1Content, link2Content, node1Content,
                         node2Content);
}