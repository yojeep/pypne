#ifndef MEDIALAXIS_H
#define MEDIALAXIS_H

void medialSurface::calc_distmap(voxel *vit, unsigned char vValue, const voxelImage &vxls, std::vector<std::vector<node>> &oldAliens) const
{
	const int i = vit->i, j = vit->j, k = vit->k;

	// 预先计算常用值
	const int nx_half = nx / 2;
	const int ny_half = ny / 2;
	const int nz_half = nz / 2;
	const int i_minus_nx_half = i - nx_half;
	const int j_minus_ny_half = j - ny_half;
	const int k_minus_nz_half = k - nz_half;

	node nalien(i, j, -nz);
	int epxMax = 2 * nx;
	double frz2 = epxMax * epxMax;
	int frz1 = 2, fry1 = 0;

	// 预先计算条件判断
	const bool is_k_gt_0 = (k > 0);
	const bool is_j_gt_0 = (j > 0);

	if (is_k_gt_0)
	{
		if (vValue != vxls(i, j, k - 1))
		{
			nalien.i = i;
			nalien.j = j;
			nalien.k = k - 1;
			frz2 = 1.0;
			frz1 = 3;
		}
		else if (is_j_gt_0)
		{
			if (vValue != vxls(i, j - 1, k))
			{
				nalien.i = i;
				nalien.j = j - 1;
				nalien.k = k;
				frz2 = 1.0;
				frz1 = 3;
			}
			else
			{
				const node &nalienOldi = oldAliens[j - 1][i];
				const int di = nalienOldi.i - i;
				const int dj = nalienOldi.j - j;
				const int dk = nalienOldi.k - k;
				const int neilienDistSqr = di * di + dj * dj + dk * dk;

				frz2 = static_cast<double>(neilienDistSqr);
				frz1 = -static_cast<int>(sqrt(neilienDistSqr)) - 1;
				fry1 = dj - 1;
				nalien = nalienOldi;
			}
		}

		const node &nalienOldi = oldAliens[j][i];
		const int di = nalienOldi.i - i;
		const int dj = nalienOldi.j - j;
		const int dk = nalienOldi.k - k;
		const int neilienDistSqr = di * di + dj * dj + dk * dk;

		if (neilienDistSqr < frz2)
		{
			nalien = nalienOldi;
			frz2 = static_cast<double>(neilienDistSqr);
			frz1 = dk - 1;
			if (j == 0)
				fry1 = -static_cast<int>(sqrt(frz2)) - 1;
		}
	}
	else if (is_j_gt_0)
	{
		if (vValue != vxls(i, j - 1, k))
		{
			nalien.i = i;
			nalien.j = j - 1;
			nalien.k = k;
			frz2 = 1.0;
			frz1 = 3;
		}
		else
		{
			const node &nalienOldi = oldAliens[j - 1][i];
			const int di = nalienOldi.i - i;
			const int dj = nalienOldi.j - j;
			const int dk = nalienOldi.k - k;
			const int neilienDistSqr = di * di + dj * dj + dk * dk;

			frz2 = static_cast<double>(neilienDistSqr);
			frz1 = 0;
			fry1 = dj - 1;
			nalien = nalienOldi;
		}
	}
	else
	{
		if (isInside(nextSegg(i, j, k).start))
		{
			epxMax = nextSegg(i, j, k).start - i;
			if (isInside(segg(i, j, k).start - 1) && i - (segg(i, j, k).start - 1) < epxMax)
			{
				epxMax = i - (segg(i, j, k).start - 1);
				nalien.i = segg(i, j, k).start - 1;
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
			nalien.i = segg(i, j, k).start - 1;
			nalien.j = j;
			nalien.k = k;
		}
		else if (isInside(i))
		{
			epxMax = std::min({nx, ny, nz}) + 1;
		}
		else
		{
			std::cout << "\n\n Error: outside voxel \n\n";
		}
		frz2 = static_cast<double>(epxMax * epxMax);
		frz1 = -epxMax;
		fry1 = -epxMax;
	}

	if (epxMax <= 0)
	{
		std::cout << i << " " << j << " " << k << "    " << segg(i, j, k).start << " " << nextSegg(i, j, k).start << std::endl;
	}

	// 优化循环计算
	const int c_start = std::max(frz1 - 1, -k);
	const int c_end = std::min(static_cast<int>(sqrt(frz2)) + 1, nz - k - 1);

	for (int c = c_start; c <= c_end; ++c)
	{
		const double c_sqr = c * c;
		if (c_sqr > frz2)
			continue;

		const double remaining_dist = frz2 - c_sqr;
		const double sqrt_remaining = sqrt(remaining_dist);

		const int b_start = std::max(std::max(static_cast<int>(-sqrt_remaining), fry1) - 1, -j);
		const int b_end = std::min(static_cast<int>(sqrt_remaining + 1.001), ny - j - 1);

		for (int b = b_start; b <= b_end; ++b)
		{
			if (vValue != vxls(i, j + b, k + c))
			{
				const int b_c_sqr = b * b + c_sqr;
				if (b_c_sqr < frz2)
				{
					frz2 = b_c_sqr;
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
					const int a = s.start - 1 - i;
					const int dist_sqr = a * a + b * b + c_sqr;
					if (dist_sqr < frz2)
					{
						frz2 = dist_sqr;
						nalien.i = i + a;
						nalien.j = j + b;
						nalien.k = k + c;
					}
				}

				const segment *next_s = &s + 1;
				if (next_s->start < nx)
				{
					const int a = next_s->start - i;
					const int dist_sqr = a * a + b * b + c_sqr;
					if (dist_sqr < frz2)
					{
						frz2 = dist_sqr;
						nalien.i = i + a;
						nalien.j = j + b;
						nalien.k = k + c;
					}
				}
			}
		}
	}

	// 边界处理
	if (!isInside(nalien.i, nalien.j, nalien.k))
	{
		nalien.i = (i_minus_nx_half < 0) ? -nx / 4 - 1 : nx * 5 / 4 + 1;
		nalien.j = (j_minus_ny_half < 0) ? -ny / 4 - 1 : ny * 5 / 4 + 1;
		nalien.k = (k_minus_nz_half < 0) ? -nz / 4 - 1 : nz * 5 / 4 + 1;

		const int dx = nalien.i - i;
		const int dy = nalien.j - j;
		const int dz = nalien.k - k;
		vit->R = sqrt(dx * dx + dy * dy + dz * dz) - 0.5;
	}
	else
	{
		const int dx = abs(nalien.i - i);
		const int dy = abs(nalien.j - j);
		const int dz = abs(nalien.k - k);

		double limit = sqrt(dx * dx + dy * dy + dz * dz) - 0.5;

		// 边界限制计算
		double iSqr = std::min(j + 2, ny - j + 1);
		if (iSqr < limit)
			limit = std::max((1.0 - _clipROutyz) * limit + _clipROutyz * iSqr, 0.01);

		iSqr = std::min(k + 2, nz - k + 1);
		if (iSqr < limit)
			limit = std::max((1.0 - _clipROutyz) * limit + _clipROutyz * iSqr, 0.01);

		iSqr = std::min(i + 2, nx - i + 1);
		if (iSqr < limit)
			limit = std::max((1.0 - _clipROutx) * limit + _clipROutx * iSqr, 0.1);

		vit->R = limit;

		// 错误检查（内联处理）
		if (frz2 <= 0)
		{
			std::cout << "WTF frz2 = " << frz2 << std::endl;
		}

		if (nalien.i < -2000 || limit > 16000000)
		{
			std::cout << "Error i = " << nalien.i << std::endl;
			std::cout << "frz2 " << frz2 << std::endl;
			std::cout << "frz1 " << frz1 << std::endl;
			std::cout << "i " << i << "  j " << j << "  k " << k << std::endl;
			std::cout << "oldAliens[j][i]. i " << oldAliens[j][i].i << "  j " << oldAliens[j][i].j << "  k " << oldAliens[j][i].k << std::endl;
			exit(0);
		}
	}

	oldAliens[j][i] = nalien;
}

voxelImage segToVxlMesh(const medialSurface &ref)
{ /// converts segments back to voxelImage
	voxelImage vxls(ref.nx, ref.ny, ref.nz, 255);
	for (int iz = 0; iz < ref.nz; ++iz)
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

	// OMPragma("omp parallel reduction(+:rBalls)")
	// std::vector<std::vector<node>> oldAliens(ny + 1, std::vector<node>(nx));
	// for (int j = 0; j < ny + 1; ++j)
	// 	for (int i = 0; i < nx; ++i)
	// 	{
	// 		oldAliens[j][i].i = i;
	// 		oldAliens[j][i].j = j;
	// 		oldAliens[j][i].k = -nz / 2 - 1;
	// 	}
	const int k = -nz / 2 - 1;
	std::vector<std::vector<node>> oldAliens;
	oldAliens.reserve(ny + 1); // 预分配外层
	for (int j = 0; j < ny + 1; ++j)
	{
		oldAliens.emplace_back();	  // 添加一个空 vector
		oldAliens.back().reserve(nx); // 预分配内层
		for (int i = 0; i < nx; ++i)
		{
			oldAliens.back().emplace_back(i, j, k); // 直接构造 node
		}
	}
	size_t nvxls10th = max(10 * int(vxlSpace.size() / 200), 1);
	// const voxel *const vnd = &*vxlSpace.end();
	const voxel *vnd = vxlSpace.data() + vxlSpace.size();
	for (voxel *vit = vxlSpace.data(); vit < vnd; ++vit)
	{
		calc_distmap(vit, 0, vxls, oldAliens);
		if (size_t(vit) % nvxls10th == 0)
		{
			(cout << "\r  distance map / sphere radius = " << vit->R).flush();
		}
		rBalls += vit->R;
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

#endif
