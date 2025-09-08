
#include "inputData.h"
#include "medialSurf.h"
#include "typses.h"

// #include "medialRadius.cpp"

medialSurface::medialSurface(inputDataNE &cfg) //, double vmvLimRelF, double crossAreaf
	: cg_(cfg), segs_(cfg.segs_), ToBeAssigned(0)
{
	setDefaults(-5.); // set _minRp to negative value if not provided by user, to be re-assigned in calc_distmaps(),  to be synced with setDefaults()

	nx = cfg.nx;
	ny = cfg.ny;
	nz = cfg.nz;

	size_t nvoxls = 0; // local copy for omp
	OMPragma("omp parallel for reduction(+:nvoxls)") for (short k = 0; k < nz; ++k) for (short j = 0; j < ny; ++j)
	{
		const segments &s = cg_.segs_[k][j];
		for (short ix = 0; ix < s.cnt; ++ix)
			if (s.s[ix].value == 0)
				nvoxls += s.s[ix + 1].start - s.s[ix].start;
	}
	nVxls = nvoxls;
	invalidSeg.start = -10000;
	invalidSeg.value = 255;
}

void medialSurface::setDefaults(double avgR)
{
	/// minRPore/Rnoise is a different keyword, of its own, not part of medialSurfaceSettings.    It is supposed to be the only adjustable parameter for the regular users.       The default value of minRPore is 1.75.

	/// medialSurfaceSettings are for advanced users and its defaults are chosen (if not provided by the user) based on the minRPore value.

	/// To extract a network that has a lower network coordination number, you can decrease the values of lenNf (e.g. 0.4), and vmvRadRelNf (e.g. 1.05) and increase the values of nRSmoothing (e.g. 9), RCorsf (e.g. 0.2) and RCors (e.g. 2.5).  You can also consider increasing minRPore (also named Rnoise) keyword to let say 2..     Every change you make you need to check that the network produces reasonable results as these are sensitive parameters and do not behave linearly. We should not be woried about the high coordination number as long as the network predicts the physical properties correctly, but not everybody agrees with me here!

	_minRp = min(1.25, avgR * 0.25) + 0.5;
	if (cg_.giv("Rnoise" + _s(0), _minRp) || cg_.giv("minRPore", _minRp) || cg_.giv("Rnoise", _minRp))
		cout << " minimum pore radius: " << _minRp << endl;
	else
		cout << " keyword \"minRPore\" not found, default value (" << abs(_minRp) << ") will be used" << endl;

	_clipROutx = 0.05;
	_clipROutyz = 0.98;
	_midRf = 0.7;
	_MSNoise = 1. * abs(_minRp) + 1.;
	_lenNf = 0.6;
	_vmvRadRelNf = 1.1;
	_nRSmoothing = 3;
	_RCorsnf = 0.15;
	_RCorsn = abs(_minRp);

	if (cg_.nBP6 == 6)
		_clipROutyz = _clipROutx;

	std::istringstream keywrdData;
	if (cg_.giv("medialSurfaceSettings" + _s(0), keywrdData) || cg_.giv("medialSurfaceSettings", keywrdData))
	{
		keywrdData >> _clipROutx >> _clipROutyz >> _midRf >> _MSNoise >> _lenNf >> _vmvRadRelNf >> _nRSmoothing >> _RCorsnf >> _RCorsn;
	}

	if (_minRp < 0.)
		cout << " Default setting, will be updated after distance map computation:\n";
	cout << "  minRPore     : " << abs(_minRp) << ";\n";
	cout << "  medialSurfaceSettings: " << _clipROutx << "  " << _clipROutyz << "  " << _midRf << "  " << _MSNoise << "  " << _lenNf << "  " << _vmvRadRelNf << "  " << _nRSmoothing << "  " << _RCorsnf << "  " << _RCorsn << endl;

	cout << "  medialSurfaceSettings:\n"
		 << "   clipROutx     : " << _clipROutx << "\n"
		 << "   clipROutyz    : " << _clipROutyz << "\n"
		 << "   midRFrac      : " << _midRf << "\n"
		 << "   RMedSurfNoise : " << _MSNoise << "\n"
		 << "   lenNf         : " << _lenNf << "\n"
		 << "   vmvRadRelNf   : " << _vmvRadRelNf << "\n"
		 << "   nRSmoothing   : " << _nRSmoothing << "\n"
		 << "   RCorsf  : " << _RCorsnf << "\n"
		 << "   RCors   : " << _RCorsn << "\n"
		 << endl;
}

void medialSurface::buildvoxelspace()
{ ///  Build voxelspace -- memory for void/active voxels
	cout << "\nProcessing " << cg_._rockTypes[0].name << " voxels:" << endl;
	cout << " Creating " << nVxls << " voxels with index: " << int(0);
	cout.flush();

	vxlSpace.resize(nVxls);

	std::vector<voxel>::iterator p = vxlSpace.begin();
	const auto vxlBegin = vxlSpace.begin();
	iZ.resize(nz);
	for (int iz = 0; iz < nz; ++iz)
	{
		for (int iy = 0; iy < ny; ++iy)
		{
			const segments &s = segs_[iz][iy];
			for (int ix = 0; ix < s.cnt; ++ix)
			{
				if (s.s[ix].value == 0)
					for (int i = s.s[ix].start; i < s.s[ix + 1].start; ++i)
					{
						p->i = i;
						p->j = iy;
						p->k = iz;
						size_t currentIndex = p - vxlBegin;
						iZ[iz].push_back(currentIndex);
						++p;
					}
			}
		}
	}

	if (nVxls != size_t(p - vxlSpace.begin()))
		cout << "\n Error created " << size_t(p - vxlSpace.begin()) << " voxels " << endl;

	cout << endl;

	/// Link voxels to segments
	p = vxlSpace.begin();
	for (int iz = 0; iz < nz; ++iz)
		for (int iy = 0; iy < ny; ++iy)
		{
			segments &s = segs_[iz][iy];
			for (int ix = 0; ix < s.cnt; ++ix)
				if (s.s[ix].value == 0)
				{
					s.s[ix].segV = &*p;
					p += s.s[ix + 1].start - s.s[ix].start;
				}
		}
}

void medialSurface::paradox_pre_removeincludedballI() // to remove the included maximal bals
{													  /// Remove maximal-balls, leave one in each adjacent voxesl. This saves time when sorting in paradoxremoveincludedballI()
	if (!nVxls)
	{
		return;
	}

	cout << " pre-remove included balls: out of " << vxlSpace.size();
	cout.flush();

	int ndel = 0;

	for (int kk = 0; kk < nz; kk += 2)
	{
		for (int jj = 0; jj < ny; jj += 2)
		{
			const segments &s = segs_[kk][jj];
			for (int ix = 0; ix < s.cnt; ++ix)
			{
				if (s.s[ix].value == 0)
				{
					for (int ii = s.s[ix].start; ii < s.s[ix + 1].start; ii += 2)
					{
						voxel *smallers[8] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
						int counter = -1;
						float maxRRR = 0;
						voxel *maxRPV = nullptr;

						for (int c = 0; c < 2; ++c)
							for (int b = 0; b < 2; ++b)
								for (int a = 0; a < 2; ++a)
								{
									voxel *vi = vxl(ii + a, jj + b, kk + c);
									if (vi != nullptr && vi->ball == &ToBeAssigned)
									{
										if (vi->R > maxRRR)
										{
											if (maxRPV)
												smallers[++counter] = maxRPV;
											maxRRR = vi->R;
											maxRPV = vi;
										}
										else
										{
											smallers[++counter] = vi;
										}
									}
								}
						++counter;
						ndel += counter;
						while (counter > 0)
							if (smallers[--counter])
							{
								smallers[counter]->ball = nullptr;
							}
					}
				}
			}
		}
	}

	nBalls -= ndel;
	cout << ",   removed = " << ndel << " remained = " << nBalls << endl;
}

void medialSurface::paradoxremoveincludedballI()
{ /// Remove included balls.  What remains are called maximal spheres
	if (!nVxls)
	{
		return;
	}

	cout << " sorting... ";
	std::vector<voxel *> tvs;
	tvs.reserve(nBalls); // allocate memory for void voxels
	{
		std::vector<voxel>::iterator ti = vxlSpace.begin() - 1;
		std::vector<voxel>::iterator tend = vxlSpace.end();
		while (++ti < tend)
			if (ti->ball)
				tvs.push_back(&(*(ti)));
	}
	cout << tvs.size() << " balls" << endl;
	sort(tvs.begin(), tvs.end(), metaballcomparer());

	cout << " remove included balls:";
	cout.flush();

	int ndel = 0;
	auto vpp = tvs.begin(), end = tvs.end();
	while (vpp < end)
	{

		voxel *vi = *vpp;
		if (!vi->ball)
		{
			++vpp;
			continue;
		}

		const int x = vi->i;
		const int y = vi->j;
		const int z = vi->k;
		const float ri = vi->R;
		const float ripinc = ri + 0.55; //.+RPreDelete
		const float mbmbDist = _RCorsnf * ri + _RCorsn;

		int ex, ey, ez;
		ex = ripinc;
		for (int a = -ex; a <= ex; ++a)
		{
			float arg_ey = ripinc * ripinc - a * a;
			if (arg_ey < 0)
				continue;
			ey = std::sqrt(arg_ey);
			for (int b = -ey; b <= ey; ++b)
			{
				float arg_ez = ripinc * ripinc - a * a - b * b;
				if (arg_ez < 0)
					continue;
				ez = sqrt(arg_ez); // sqrts(r2i)+1-a-b;
				for (int c = -ez; c <= ez; ++c)
				{
					voxel *vj = vxl(x + a, y + b, z + c);
					if ((vj != nullptr) && (vj->ball) && (vj != vi))
					{
						const float rj = vj->R;
						if (rj <= ri)
						{
							float D = sqrtf(a * a + b * b + c * c);

							if (D < mbmbDist || (D + rj < ripinc + _MSNoise))
							{
								vj->ball = nullptr;
								++ndel;
							}
						}
					}
				}
			}
		}

		++vpp;
		if ((vpp - tvs.begin()) % 10000 == 0)
			cout << "\r  remove = " << ndel;
	}
	cout << "\r  removed = " << ndel << " remained = " << tvs.size() - ndel << " balls" << endl;
	nBalls -= ndel;
	return;
}

void medialSurface::moveUphill(medialBall *b_i) // const
{												/// Refines the maximal-sphere location and radius,

	const voxel *vi = vxl(b_i->fi, b_i->fj, b_i->fk);
	dbl3 disp(0., 0., 0.);
	{
		const voxel *vjm = vxl(vi->i - 1, vi->j, vi->k);
		const voxel *vjp = vxl(vi->i + 1, vi->j, vi->k);
		if (vjm && vjp)
		{
			float gp = vjp->R - vi->R;
			float gm = vi->R - vjm->R;
			if (abs(gp - gm) > 0.01)
				disp.x = max(-0.49, min(0.49, -0.5 * (gp + gm) / (gp - gm)));
		}
	}
	{
		const voxel *vjm = vxl(vi->i, vi->j - 1, vi->k);
		const voxel *vjp = vxl(vi->i, vi->j + 1, vi->k);
		if (vjm && vjp)
		{
			float gp = vjp->R - vi->R;
			float gm = vi->R - vjm->R;
			if (abs(gp - gm) > 0.01)
				disp.y = max(-0.49, min(0.49, -0.5 * (gp + gm) / (gp - gm)));
		}
	}
	{
		const voxel *vjm = vxl(vi->i, vi->j, vi->k - 1);
		const voxel *vjp = vxl(vi->i, vi->j, vi->k + 1);
		if (vjm && vjp)
		{
			float gp = vjp->R - vi->R;
			float gm = vi->R - vjm->R;
			if (abs(gp - gm) > 0.01)
				disp.z = max(-0.49, min(0.49, -0.5 * (gp + gm) / (gp - gm)));
		}
	}
	if (b_i != b_i->boss)
	{
		dbl3 BosKidVec = *b_i - *(b_i->boss);
		disp -= 0.95 * ((BosKidVec & disp) / (magSqr(BosKidVec) + 1e-12)) * BosKidVec;
	}
	b_i->fi = vi->i - _mp5 + disp.x;
	b_i->fj = vi->j - _mp5 + disp.y;
	b_i->fk = vi->k - _mp5 + disp.z;
	b_i->R = vi->R + 0.95 * mag(disp);
}

void medialSurface::moveUphillp1(medialBall *bi) // const
{												 /// Refines the maximal-sphere location, by moving it uphil the gradient of the distance-map, potentially relocates to new voxels

	const voxel *vi = vxl(bi->fi, bi->fj, bi->fk);
	dbl3 disp(0., 0., 0.), grad(0., 0., 0.);

	{
		const voxel *vjm = vxl(vi->i - 1, vi->j, vi->k);
		const voxel *vjp = vxl(vi->i + 1, vi->j, vi->k);
		if (vjm && vjp)
		{
			float gp = vjp->R - vi->R;
			float gm = vi->R - vjm->R;
			grad.x = 0.5 * (gp + gm);
			if (abs(gp - gm) > 0.01)
				disp.x = max(-0.59, min(0.59, -0.5 * (gp + gm) / (gp - gm)));
		}
	}
	{
		const voxel *vjm = vxl(vi->i, vi->j - 1, vi->k);
		const voxel *vjp = vxl(vi->i, vi->j + 1, vi->k);
		if (vjm && vjp)
		{
			float gp = vjp->R - vi->R;
			float gm = vi->R - vjm->R;
			grad.y = 0.5 * (gp + gm);
			if (abs(gp - gm) > 0.01)
				disp.y = max(-0.59, min(0.59, -0.5 * (gp + gm) / (gp - gm)));
		}
	}
	{
		const voxel *vjm = vxl(vi->i, vi->j, vi->k - 1);
		const voxel *vjp = vxl(vi->i, vi->j, vi->k + 1);
		if (vjm && vjp)
		{
			float gp = vjp->R - vi->R;
			float gm = vi->R - vjm->R;
			grad.z = 0.5 * (gp + gm);
			if (abs(gp - gm) > 0.01)
				disp.z = max(-0.59, min(0.59, -0.5 * (gp + gm) / (gp - gm)));
		}
	}
	disp += 1.4 * grad;

	if (bi != bi->boss)
	{
		dbl3 BosKidVec = *bi - *(bi->boss);
		disp -= 0.5 * ((BosKidVec & disp) / (magSqr(BosKidVec) + 1e-12)) * BosKidVec;
	}
	disp /= (0.55 * mag(disp) + 0.05);

	voxel *vxlj = vxl(bi->fi + disp[0], bi->fj + disp[1], bi->fk + disp[2]);
	if (vxlj && vi != vxlj && vxlj->R > vi->R && vxlj->ball == nullptr)
	{
		bi->fi = vxlj->i - _mp5;
		bi->fj = vxlj->j - _mp5;
		bi->fk = vxlj->k - _mp5;
		bi->R = vxlj->R;
		bi->vxl->ball = nullptr;
		bi->vxl = vxlj;
		vxlj->ball = bi;
		//++nrelocations;
	}

	// cout<<nrelocations<<" relocations "<<endl;
}

void makeFriend(medialBall *vi, medialBall *vj)
{
	if (vj->R > vi->R)
	{
		medialBall *tmp = vi;
		vi = vj;
		vj = tmp;
	}
	if ((vi->R < 1.5 * vj->R) && (!vi->isNei(vj)) && (!vi->inParents(vj)) && !vj->inParents(vi))
	{
		vi->addNei(vj);
		vj->addNei(vi);
	}
}

/*inline double cosAngleWithBossPerD(const medialBall* a, const medialBall* b)  {
	dbl3 v1=*a-*b;
	if (b==b->boss) return 1./(mag(v1));
	dbl3 v2=*b-*(b->boss);
	double dotProd=v2&v1;
	return sqrt(dotProd*dotProd/(magSqr(v1)*magSqr(v1)*magSqr(v2)));
}*/

void medialSurface::competeForParent(medialBall *vi, medialBall *vj)
{
	const double noise = _MSNoise;

	const double ri = vi->R;
	const double rj = vj->R;
	const double riSqr = ri * ri;
	const double rjSqr = rj * rj;
	const double dSqr = distSqr(vi, vj);
	const double dVal = sqrt(dSqr);

	const double wsinv = 1.0 / (riSqr + rjSqr);
	const voxel *middlevxl = vxl(wsinv * (vi->fi * rjSqr + vj->fi * riSqr),
								 wsinv * (vi->fj * rjSqr + vj->fj * riSqr),
								 wsinv * (vi->fk * rjSqr + vj->fk * riSqr));

	if (!middlevxl)
		return;

	const double minR = min(ri, rj);
	if (middlevxl->R <= minR * _midRf - 0.5)
		return;
	if (1.01 * dVal >= ri + rj + 1.0 + noise)
		return;

	// Handle boss assignment cases
	if (vj->boss == vj && vi->mastrSphere() != vj)
	{
		if (ri >= rj)
			vj->boss = vi;
		else if (vi->boss->R <= rj)
			vi->boss = vj;
		else if (ri >= rj - noise && ri * _vmvRadRelNf + noise >= rj)
			vj->boss = vi;
	}
	else if (vi->boss == vi && vj->mastrSphere() != vi)
	{
		if (rj >= ri)
			vi->boss = vj;
		else if (vj->boss->R <= ri)
			vj->boss = vi;
		else if (rj >= ri - noise && rj * _vmvRadRelNf + noise >= ri)
			vi->boss = vj;
	}

	medialBall *mvi = vi->mastrSphere();
	medialBall *mvj = vj->mastrSphere();

	if (mvi == vj || mvj == vi)
		return;

	if (mvi == mvj)
	{
		const short leveli = vi->level();
		const short levelj = vj->level();
		const short levelDiff = leveli - levelj;

		const double distViBoss = dist(vi->boss, vi);
		const double distVjBoss = dist(vj->boss, vj);
		const double distViVj = dist(vi, vj);

		if (levelDiff < -1) // leveli + 1 < levelj
		{
			if ((vj->boss->R - vj->R + 2 * noise) / (distVjBoss + 0.25) <
				(vi->R - vj->R + 2 * noise + 0.01) / (distViVj + 0.2))
				vj->boss = vi;
		}
		else if (levelDiff > 1) // leveli > levelj + 1
		{
			if ((vi->boss->R - vi->R + 2 * noise) / (distViBoss + 0.25) <
				(vj->R - vi->R + 2 * noise + 0.01) / (distViVj + 0.2))
				vi->boss = vj;
		}
		else
		{
			if (levelDiff > 0) // leveli > levelj
			{
				if ((vi->boss->R - vi->R + 2 * noise) / (distViBoss + 1.2) <
						(vj->R - vi->R + 2 * noise) / (distViVj + 1.3) &&
					!vj->inParents(vi))
					vi->boss = vj;
			}
			else if (levelDiff < 0) // leveli < levelj
			{
				if ((vj->boss->R - vj->R + 2 * noise) / (distVjBoss + 1.2) <
						(vi->R - vj->R + 2 * noise) / (distViVj + 1.3) &&
					!vi->inParents(vj))
					vj->boss = vi;
			}

			if (middlevxl->R >= 0.45 * (ri + rj) - 1.0 && dVal < (ri + rj) * 0.5 + 2.0)
				makeFriend(vi, vj);
		}
	}
	else // mvi != mvj
	{
		const double avgR = 0.5 * (mvi->R + mvj->R);
		if (distSqr(mvi, mvj) <= _lenNf * (avgR + 2 * noise) * (avgR + 2 * noise))
		{
			// Ensure mvi is the larger one
			if (mvi->R < mvj->R)
			{
				std::swap(vi, vj);
				std::swap(mvi, mvj);
			}

			if (mvj->R < _vmvRadRelNf * vj->R + noise &&
				mvj->R < _vmvRadRelNf * vi->R + noise &&
				mvj->R < _vmvRadRelNf * vi->boss->R + noise)
			{
				while (vj != vj->boss && mvj->R < _vmvRadRelNf * vj->boss->R + noise)
				{
					medialBall *pvj = vj->boss;
					vj->boss = vi;
					vi = vj;
					vj = pvj;
				}
				if (vj->boss == vj && vi->mastrSphere() != vj)
					vj->boss = vi;
			}
		}

		if (vi != vj->boss)
		{
			mvi = vi->mastrSphere();
			mvj = vj->mastrSphere();

			short leveli = vi->level();
			short levelj = vj->level();
			const double distMviMvj = dist(mvi, mvj);
			const double distAvg = distMviMvj + 0.5 * noise;

			while (leveli >= levelj)
			{
				const double viBossRatio = (vi->boss->R - vi->R + 0.55 * noise) / (dist(mvi, vi) + distAvg);
				const double vjRatio = (vj->R - vi->R + 0.5 * noise) / (dist(mvj, vi) + distAvg);
				if (viBossRatio >= vjRatio)
					break;

				medialBall *pvi = vi->boss;
				vi->boss = vj;
				vj = vi;
				vi = pvi;
				++levelj;
				--leveli;
			}

			while (levelj >= leveli)
			{
				const double vjBossRatio = (vj->boss->R - vj->R + 0.55 * noise) / (dist(mvj, vj) + distAvg);
				const double viRatio = (vi->R - vj->R + 0.5 * noise) / (dist(mvi, vj) + distAvg);
				if (vjBossRatio >= viRatio)
					break;

				medialBall *pvj = vj->boss;
				vj->boss = vi;
				vi = vj;
				vj = pvj;
				++leveli;
				--levelj;
			}

			makeFriend(vi, vj);
		}
	}
}

void medialSurface::findBoss(medialBall *vi)
{

	const float x = vi->fi, y = vi->fj, z = vi->fk;
	const float ripp = vi->R * 0.6 + 2. * _MSNoise + 2.;
	const float ex = x + ripp;
	for (float xpa = 2. * x - ex; xpa <= ex; xpa += 1.0f)
	{
		float remain_y = ripp * ripp - (xpa - x) * (xpa - x);
		if (remain_y < 0.f)
			continue;
		float ey = y + sqrt(remain_y);
		for (float ypb = 2. * y - ey; ypb <= ey; ypb += 1.0f)
		{
			float remain_z = remain_y - (y - ypb) * (y - ypb);
			if (remain_z < 0.f)
				continue;
			float ez = z + sqrt(remain_z);
			for (float zpc = 2. * z - ez; zpc <= ez; zpc += 1.0f)
			{
				voxel *vj = this->vxl(xpa, ypb, zpc);
				if ((vj != nullptr) && vj->ball && (*vi != *vj))
				{ //--------------------------------------------------------
					competeForParent(&*vi, vj->ball);
				} //--------------------------------------------------------
			}
		}
	}
}

voxelImage segToVxlMesh(const medialSurface &ref)
{ /// converts segments back to voxelImage
	voxelImage vxls(ref.nx, ref.ny, ref.nz, 255);
	OMPFor() for (int iz = 0; iz < ref.nz; ++iz)
	{
		for (int iy = 0; iy < ref.ny; ++iy)
		{
			const segments &s = ref.segs_[iz][iy];
			for (int ix = 0; ix < s.cnt; ++ix)
			{
				std::fill(&vxls(s.s[ix].start, iy, iz), &vxls(s.s[ix + 1].start, iy, iz), s.s[ix].value);
			}
		}
	}
	return vxls;
}

void medialSurface::calc_distmap(voxel &vit, unsigned char vValue, const voxelImage &vxls, std::vector<std::vector<node>> &oldAliens) const
{

	const int i = vit.i, j = vit.j, k = vit.k;

	node nalien(i, j, -nz);

	int
		epxMax = 2 * nx;

	double frz2 = epxMax * epxMax;
	int frz1 = 2, fry1 = 0;

	if (k > 0)
	{
		if (vValue != vxls(i, j, k - 1))
		{
			nalien.i = i;
			nalien.j = j;
			nalien.k = k - 1;
			frz2 = 1.;
			frz1 = 3;
			// fry1 = 0;
		}
		else
		{

			if (j > 0)
			{
				if (vValue != vxls(i, j - 1, k))
				{
					nalien.i = i;
					nalien.j = j - 1;
					nalien.k = k;
					frz2 = 1.;
					frz1 = 3;
					// fry1 = 0;
				}
				else
				{
					const node &nalienOldi = oldAliens[j - 1][i];

					int neilienDistSqr = (nalienOldi.i - i) * (nalienOldi.i - i) + (nalienOldi.j - j) * (nalienOldi.j - j) + (nalienOldi.k - k) * (nalienOldi.k - k);
					frz2 = neilienDistSqr + 0.;
					frz1 = -sqrt(neilienDistSqr) - 1;
					fry1 = nalienOldi.j - j - 1;
					nalien.i = nalienOldi.i;
					nalien.j = nalienOldi.j;
					nalien.k = nalienOldi.k;
				}
			}

			const node &nalienOldi = oldAliens[j][i];
			int neilienDistSqr = (nalienOldi.i - i) * (nalienOldi.i - i) + (nalienOldi.j - j) * (nalienOldi.j - j) + (nalienOldi.k - k) * (nalienOldi.k - k);
			if (neilienDistSqr < frz2)
			{
				nalien.i = nalienOldi.i;
				nalien.j = nalienOldi.j;
				nalien.k = nalienOldi.k;

				frz2 = neilienDistSqr + 0.;
				frz1 = nalienOldi.k - k - 1;
				if (j == 0)
					fry1 = -sqrt(frz2) - 1;
			}
		}
	}
	else if (j > 0)
	{
		if (vValue != vxls(i, j - 1, k))
		{
			nalien.i = i;
			nalien.j = j - 1;
			nalien.k = k;
			frz2 = 1.;
			frz1 = 3;
			// fry1 = 0;
		}
		else
		{
			const node &nalienOldi = oldAliens[j - 1][i];

			int neilienDistSqr = (nalienOldi.i - i) * (nalienOldi.i - i) + (nalienOldi.j - j) * (nalienOldi.j - j) + (nalienOldi.k - k) * (nalienOldi.k - k);
			frz2 = neilienDistSqr + 0.;
			frz1 = 0;
			fry1 = nalienOldi.j - j - 1;
			nalien.i = nalienOldi.i;
			nalien.j = nalienOldi.j;
			nalien.k = nalienOldi.k;
		}
	}
	else // if(	epxMax == 2*nx)//. XXXXXX WARNING
	{
		if (isInside(nextSegg(i, j, k).start))
		{
			epxMax = nextSegg(i, j, k).start - i;
			if (isInside(segg(i, j, k).start - 1) && i - (segg(i, j, k).start - 1) < epxMax)
			{
				epxMax = i - (segg(i, j, k).start - 1);
				nalien.i = (segg(i, j, k).start - 1);
				nalien.j = j;
				nalien.k = k;
			}
			else
			{
				nalien.i = nextSegg(i, j, k).start;
				nalien.j = j;
				nalien.k = k;
			}
		}
		else if (isInside(segg(i, j, k).start - 1))
		{
			epxMax = std::min(i - segg(i, j, k).start + 1, epxMax);
			nalien.i = (segg(i, j, k).start - 1);
			nalien.j = j;
			nalien.k = k;
		}
		else if (isInside(i))
		{
			epxMax = std::min(nx, std::min(ny, nz)) + 1;
		}
		else
		{
			cout << "\n\n Error: outside voxel \n\n";
		}
		frz2 = epxMax * epxMax + 0.;
		frz1 = -epxMax;
		fry1 = -epxMax;
	}

	if (epxMax <= 0)
		cout << i << " " << j << " " << k << "    " << segg(i, j, k).start << " " << (nextSegg(i, j, k)).start << " " << endl;

	for (int c = std::max(frz1 - 1, -k); c <= min(int(sqrt(frz2)) + 1, nz - k - 1); ++c)
	{
		// if( isInside(i, j, k+c))
		{
			const int blim = min(int(sqrt(frz2 - c * c) + 1.001), ny - j - 1);
			for (int b = std::max(std::max(int(-sqrt(frz2 - c * c)), fry1) - 1, -j); b <= blim; ++b)
			{

				// if(isJInside(j+b)) //
				{
					if (vValue != vxls(i, j + b, k + c))
					{
						if ((b * b + c * c) < frz2)
						{
							frz2 = b * b + c * c;
							nalien.i = i;
							nalien.j = j + b;
							nalien.k = k + c;
						}
					}
					else
					{
						const segment &s = segg(i, j + b, k + c);
						if (s.start > 0)
						{
							int a = (s.start - 1 - i);
							if ((a * a + b * b + c * c) < frz2)
							{
								frz2 = a * a + b * b + c * c;
								nalien.i = i + a;
								nalien.j = j + b;
								nalien.k = k + c;
							}
						}

						if ((&s + 1)->start < nx)
						{
							int a = ((&s + 1)->start - i);
							if ((a * a + b * b + c * c) < frz2)
							{
								frz2 = a * a + b * b + c * c;
								nalien.i = i + a;
								nalien.j = j + b;
								nalien.k = k + c;
							}
						}
					}
				}
			}
		}
	}

	if (!isInside(nalien.i, nalien.j, nalien.k))
	{
		nalien.i = (i < nx / 2) ? -nx / 4 - 1 : nx * 5 / 4 + 1;
		nalien.j = (j < ny / 2) ? -ny / 4 - 1 : ny * 5 / 4 + 1;
		nalien.k = (k < nz / 2) ? -nz / 4 - 1 : nz * 5 / 4 + 1;
		vit.R = sqrt((nalien.i - i) * (nalien.i - i) + (nalien.j - j) * (nalien.j - j) + (nalien.k - k) * (nalien.k - k)) - 0.5;
	}
	else
	{
		int dx = abs(nalien.i - i), dy = abs(nalien.j - j), dz = abs(nalien.k - k);

		double limit = sqrt(dx * dx + dy * dy + dz * dz) - 0.5;
		double iSqr = min((j + 2), (ny - j + 1));
		if (iSqr < limit)
			limit = max((1. - _clipROutyz) * limit + _clipROutyz * iSqr, 0.01);
		iSqr = min((k + 2), (nz - k + 1));
		if (iSqr < limit)
			limit = max((1. - _clipROutyz) * limit + _clipROutyz * iSqr, 0.01);
		iSqr = min((i + 2), (nx - i + 1));
		if (iSqr < limit)
			limit = max((1. - _clipROutx) * limit + _clipROutx * iSqr, 0.1);
		vit.R = limit;

		if (frz2 <= 0)
			cout << "WTF frz2 = " << frz2 << endl;
		if (nalien.i < -2000 || limit > 16000000)
		{
			cout << "Error i = " << nalien.i << endl;
			cout << "frz2 " << frz2 << endl;
			cout << "frz1 " << frz1 << endl;
			cout << "i " << i << "  j " << j << "  k " << k << endl;
			cout << "oldAliens[j][i]. i " << oldAliens[j][i].i << "  j " << oldAliens[j][i].j << "  k " << oldAliens[j][i].k << endl;
			exit(0);
		}
	}

	// OMPragma("omp critical")
	oldAliens[j][i] = nalien;
}

void medialSurface::calc_distmaps() // search  MBs at each voxel
{

	cout << " computing distance map for index " << int(0);
	cout.flush();

	if (!nVxls)
	{
		cout << " no voxels no balls,\n"
			 << endl;
		return;
	}

	voxelImage vxls = segToVxlMesh(*this);
	double rBalls = 0.;
	const int k = -nz / 2 - 1;
	thread_local std::vector<std::vector<node>> oldAliens;
	size_t nvxls10th = max(10 * int(vxlSpace.size() / 200), 1);
	// 使用iZ进行z方向并行处理
	OMPragma("omp parallel for reduction(+:rBalls) schedule(dynamic)") for (int iz = 0; iz < nz; ++iz)
	{
		// 使用静态局部变量确保每个线程只初始化一次
		static thread_local bool initialized = false;
		// cout << initialized << endl;
		if (!initialized)
		{
			oldAliens.reserve(ny + 1);
			for (int j = 0; j < ny + 1; ++j)
			{
				oldAliens.emplace_back();
				oldAliens.back().reserve(nx);
				for (int i = 0; i < nx; ++i)
				{
					oldAliens.back().emplace_back(i, j, k);
				}
			}
			initialized = true;
		}

		for (size_t idx : iZ[iz])
		{
			voxel &vit = vxlSpace[idx];
			calc_distmap(vit, 0, vxls, oldAliens);

			// 控制进度输出频率并避免线程竞争
			if (idx % nvxls10th == 0)
			{
				OMPragma("omp critical")
				{
					(cout << "\r  distance map / sphere radius = " << vit.R).flush();
				}
			}
			rBalls += vit.R;
		}
		for (size_t idx : iZ[iz])
		{
			voxel &vit = vxlSpace[idx];
			const int i = vit.i, j = vit.j;
			node &oldAliens_ = oldAliens[j][i];
			oldAliens_.i = i;
			oldAliens_.j = j;
			oldAliens_.k = k;
		}
	}
	cout << "\n  average distance map = " << rBalls / nVxls << endl;

	if (_minRp < 0.)
	{
		setDefaults(rBalls / nVxls);
	}

	return;
}

void medialSurface::smoothRadius()
{

	(cout << " smoothing R  ").flush();

	std::vector<float> delRrr(vxlSpace.size(), 0.0f);
	(cout << "*").flush();
	OMPFor() for (int k = 0; k < nz; ++k)
	{
		for (int j = 0; j < ny; ++j)
		{
			const segments &s = cg_.segs_[k][j];
			for (int ix = 0; ix < s.cnt; ++ix)
				if (s.s[ix].value == 0)
				{
					segment &seg = s.s[ix];
					for (int i = seg.start; i < s.s[ix + 1].start; ++i)
					{
						double sumR = 0.;
						int counter = 0;
						for (int kk = max(k - 1, 0); kk < min(k + 2, nz); ++kk)
							for (int jj = max(j - 1, 0); jj < min(j + 2, ny); ++jj)
							{
								int ii = max(i - 1, 0);
								const segment *segbc = cg_.segptr(ii, jj, kk);
								if (segbc->value != 0 && (segbc + 1)->value == 0)
								{
									++segbc;
									ii = segbc->start;
								}
								if (segbc->value == 0)
								{
									int ii2 = min((segbc + 1)->start, i + 2);
									voxel *vxlj = segbc->segV + (ii - segbc->start);
									for (; ii < ii2; ++ii)
									{
										sumR += vxlj->R;
										++vxlj;
										counter += 1;
									}
								}
							}

						delRrr[seg.segV + (i - seg.start) - (&vxlSpace[0])] = 4. * sumR / (3 * counter + 27) - seg.segV[i - seg.start].R;
					}
				}
		}
	}
	(cout << "*").flush();

	OMPFor() for (int k = 0; k < nz; ++k)
	{
		for (int j = 0; j < ny; ++j)
		{
			const segments &s = cg_.segs_[k][j];
			for (int ix = 0; ix < s.cnt; ++ix)
				if (s.s[ix].value == 0)
				{
					segment &seg = s.s[ix];
					for (int i = seg.start; i < s.s[ix + 1].start; ++i)
					{
						double sumDelR = 0.;
						int counter = 0;
						for (int kk = max(k - 1, 0); kk < min(k + 2, nz); ++kk)
							for (int jj = max(j - 1, 0); jj < min(j + 2, ny); ++jj)
							{
								int ii = max(i - 1, 0);
								const segment *segbc = cg_.segptr(ii, jj, kk);
								if (segbc->value != 0 && (segbc + 1)->value == 0)
								{
									++segbc;
									ii = segbc->start;
								}
								if (segbc->value == 0)
								{
									int ii2 = min((segbc + 1)->start, i + 2);
									voxel *vxlj = segbc->segV + (ii - segbc->start);
									for (; ii < ii2; ++ii)
									{
										sumDelR += delRrr[vxlj - (&vxlSpace[0])];
										++vxlj;
										++counter;
									}
								}
							}

						seg.segV[i - seg.start].R += min(max(0.02 * (delRrr[seg.segV + (i - seg.start) - (&vxlSpace[0])] - 0.99 * 2. * sumDelR / (1 * counter + 27)), -0.005), 0.01);
					}
				}
		}
	}
	(cout << "*").flush();

	{ /// Finally, report max distance map, to confirm that distance map is not changed too much
		float maxrrr = 0;
		OMPragma("omp parallel for reduction(max:maxrrr)") for (auto ti = vxlSpace.begin(); ti < vxlSpace.end(); ++ti) maxrrr = max(maxrrr, ti->R);
		cout << " maxrrr " << maxrrr << endl;
	}
}

void medialSurface::createBallsAndHierarchy()
{ /// Create distance map, maximal-spheres, and their hirarchy (medial-surface connectivity)

	//	mediaAxes medAxis(*this, minRP, clipOutSideBallFraction, clipOutSideBallFraction*0.5+0.);

	buildvoxelspace();

	calc_distmaps();

	for (int i = 0; i < _nRSmoothing; ++i)
		smoothRadius();

	nBalls = 0;
	double rBalls = 0.;
	std::vector<voxel>::iterator vit = vxlSpace.begin() - 1;
	const std::vector<voxel>::iterator vend = vxlSpace.end();
	while (++vit < vend)
	{
		if (vit->R >= _minRp)
		{
			vit->ball = &(ToBeAssigned);
			++nBalls;
			rBalls += vit->R;
		}
		else
		{
			vit->ball = nullptr;
		}
	}
	cout << "\n  number of potential maximal spheres: " << nBalls << ",  average radius = " << rBalls / nBalls << endl;

	paradox_pre_removeincludedballI();

	paradoxremoveincludedballI();

	cout << " collecting maximal balls out of " << nBalls << endl;

	std::vector<voxel *> tvs;
	tvs.reserve(nBalls);
	{
		std::vector<voxel>::iterator ti = vxlSpace.begin() - 1;
		std::vector<voxel>::iterator tend = vxlSpace.end();
		while (++ti < tend)
			if (ti->ball)
			{
				if ((ti->R) >= _minRp)
					tvs.push_back(&(*(ti)));
				else
					cout << "  sdsd ";
			}

		cout << " sorting " << int(tvs.size()) << " maximal balls" << endl;
		sort(tvs.begin(), tvs.end(), metaballcomparer());
	}

	ballSpace.reserve(nBalls);
	{
		std::vector<voxel *>::iterator ti = tvs.begin() - 1;
		std::vector<voxel *>::iterator tend = tvs.end();
		while (++ti < tend)
		{
			ballSpace.emplace_back(*ti, 0);
			(*ti)->ball = &*(ballSpace.rbegin());
		}
	}

	const std::vector<medialBall>::iterator voxend = ballSpace.end();

	{
		std::vector<medialBall>::iterator vi = ballSpace.begin() - 1;
		while (++vi != voxend)
			moveUphill(&*vi);
	}
	{
		std::vector<medialBall>::iterator vi = ballSpace.begin() - 1;
		while (++vi != voxend)
			moveUphillp1(&*vi);
	}

	{
		std::vector<medialBall>::iterator vi = ballSpace.begin() - 1;
		while (++vi != voxend)
			moveUphill(&*vi);
	}

	cout << " creating ball hierarchy:";
	cout.flush();
	const std::vector<medialBall>::iterator voxp = ballSpace.begin();
	{
		std::vector<medialBall>::iterator vi = ballSpace.begin();
		while (vi != voxend)
		{
			findBoss(&*vi);
			if ((vi - voxp) % 100000 == 0)
			{
				cout << "\r   ball: " << int(vi - voxp);
				cout.flush();
			}
			++vi;
		}
		cout << "\r   ball: " << int(vi - voxp) << endl;
	}
}