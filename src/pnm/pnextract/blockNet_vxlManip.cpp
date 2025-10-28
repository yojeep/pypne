#ifndef voxelImageManip_H
#define voxelImageManip_H

#include "blockNet.h"
// #include "global.h"
// class mapComparer  {  public:
// bool operator() (pair<const int,short>& i1, pair<const int,short> i2) {return i1.second<i2.second;}  };

using namespace std; // std::pair, vector map

size_t growPores_X2(voxelField<int> &VElems, int bgn, int lst, int porValue)
{
	auto &pool = GlobalThreadPool::get();
	size_t nChanges(0);
	const int nz = VElems.nz(), ny = VElems.ny(), nx = VElems.nx();
	{
		const voxelField<int> voxls = VElems;
		const int i_start = 1, i_end = nx - 1;
		const int j_start = 1, j_end = ny - 1;
		const int k_start = 1, k_end = nz - 1;
		const int i_size = i_end - i_start;
		const int j_size = j_end - j_start;
		const int k_size = k_end - k_start;
		const size_t total_iterations = static_cast<size_t>(k_size * j_size * i_size);
		auto nChanges_futures = pool.submit_blocks(
			0, total_iterations,
			[&](const size_t start_idx, const size_t end_idx) -> size_t
			{
				size_t local_nChanges = 0;
				for (size_t idx = start_idx; idx < end_idx; ++idx)
				{
					size_t temp = idx;
					int i = i_start + (temp % i_size);
					temp /= i_size;
					int j = j_start + (temp % j_size);
					temp /= j_size;
					int k = k_start + temp;
					const int *pijk = &voxls(i, j, k);
					if (*pijk == porValue)
					{
						if (bgn <= voxls.v_i(1, pijk) && voxls.v_i(1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_i(1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_i(-1, pijk) && voxls.v_i(-1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_i(-1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_j(1, pijk) && voxls.v_j(1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_j(1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_j(-1, pijk) && voxls.v_j(-1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_j(-1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_k(1, pijk) && voxls.v_k(1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_k(1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_k(-1, pijk) && voxls.v_k(-1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_k(-1, pijk);
							++local_nChanges;
						}
					}
				}
				return local_nChanges;
			});
		for (auto &future : nChanges_futures)
		{
			nChanges += future.get();
		}
	}

	cout << "  ngrowX3:" << nChanges << ",";
	nChanges = 0;
	{
		const voxelField<int> voxls = VElems;
		const int i_start = 1, i_end = nx - 1;
		const int j_start = 1, j_end = ny - 1;
		const int k_start = 1, k_end = nz - 1;
		const int i_size = i_end - i_start;
		const int j_size = j_end - j_start;
		const int k_size = k_end - k_start;
		const size_t total_iterations = static_cast<size_t>(k_size * j_size * i_size);
		auto nChanges_futures = pool.submit_blocks(
			0, total_iterations,
			[&](const size_t start_idx, const size_t end_idx) -> size_t
			{
				size_t local_nChanges = 0;
				for (size_t idx = start_idx; idx < end_idx; ++idx)
				{
					size_t temp = idx;
					int i = i_start + (temp % i_size);
					temp /= i_size;
					int j = j_start + (temp % j_size);
					temp /= j_size;
					int k = k_start + temp;
					const int *pijk = &voxls(i, j, k);
					if (*pijk == porValue)
					{
						if (bgn <= voxls.v_i(1, pijk) && voxls.v_i(1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_i(1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_i(-1, pijk) && voxls.v_i(-1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_i(-1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_j(1, pijk) && voxls.v_j(1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_j(1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_j(-1, pijk) && voxls.v_j(-1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_j(-1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_k(1, pijk) && voxls.v_k(1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_k(1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_k(-1, pijk) && voxls.v_k(-1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_k(-1, pijk);
							++local_nChanges;
						}
					}
				}
				return local_nChanges;
			});
		for (auto &future : nChanges_futures)
		{
			nChanges += future.get();
		}
		cout << nChanges << ",";
	}
	nChanges = 0;
	{
		const voxelField<int> voxls = VElems;
		// 内部区域
		const int i_start = 1, i_end = nx - 1; // i in [1, nx-2]
		const int j_start = 1, j_end = ny - 1; // j in [1, ny-2]
		const int k_start = 1, k_end = nz - 1; // k in [1, nz-2]
		const int i_size = i_end - i_start;
		const int j_size = j_end - j_start;
		const int k_size = k_end - k_start;
		const size_t total_iterations = static_cast<size_t>(k_size * j_size * i_size);
		auto nChanges_futures = pool.submit_blocks(
			0, total_iterations,
			[&](const size_t start_idx, const size_t end_idx) -> size_t
			{
				size_t local_nChanges = 0;
				// 我们要逆序处理：从 idx = total-1 到 0
				// 但线程池给的是 [start_idx, end_idx)
				// 所以我们让每个线程处理的块是“逻辑逆序”的
				// 方法：将线程池的正序 idx 映射为逆序的 (i,j,k)
				for (size_t idx = start_idx; idx < end_idx; ++idx)
				{
					// idx ∈ [0, total) 正序
					// 我们想要逆序：从 total-1 开始
					size_t rev_idx = total_iterations - 1 - idx; // 映射到逆序位置
					// 从 rev_idx 还原 (k, j, i)
					size_t temp = rev_idx;
					int i = i_start + (temp % i_size);
					temp /= i_size;
					int j = j_start + (temp % j_size);
					temp /= j_size;
					int k = k_start + temp;
					const int *pijk = &voxls(i, j, k);
					if (*pijk == porValue)
					{
						if (bgn <= voxls.v_i(1, pijk) && voxls.v_i(1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_i(1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_i(-1, pijk) && voxls.v_i(-1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_i(-1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_j(1, pijk) && voxls.v_j(1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_j(1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_j(-1, pijk) && voxls.v_j(-1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_j(-1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_k(1, pijk) && voxls.v_k(1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_k(1, pijk);
							++local_nChanges;
						}
						else if (bgn <= voxls.v_k(-1, pijk) && voxls.v_k(-1, pijk) <= lst)
						{
							VElems(i, j, k) = voxls.v_k(-1, pijk);
							++local_nChanges;
						}
					}
				}
				return local_nChanges;
			});
		for (auto &future : nChanges_futures)
		{
			nChanges += future.get();
		}
	}
	cout << "  ngrowX2:" << nChanges << "  ";
	return nChanges;
}

void growPores(voxelField<int> &VElems, int bgn, int lst, int porValue)
{

	auto &pool = GlobalThreadPool::get();
	const voxelField<int> voxls = VElems;
	size_t nChanges(0);
	const int nz = VElems.nz(), ny = VElems.ny(), nx = VElems.nx();
	const int i_start = 1, i_end = nx - 1;
	const int j_start = 1, j_end = ny - 1;
	const int k_start = 1, k_end = nz - 1;
	const int i_size = i_end - i_start;
	const int j_size = j_end - j_start;
	const int k_size = k_end - k_start;
	const size_t total_iterations = static_cast<size_t>(k_size * j_size * i_size);
	auto nChanges_futures = pool.submit_blocks(
		0, total_iterations,
		[&](const size_t start_idx, const size_t end_idx) -> size_t
		{
			size_t local_nChanges = 0;
			for (size_t idx = start_idx; idx < end_idx; ++idx)
			{
				size_t temp = idx;
				int i = i_start + (temp % i_size);
				temp /= i_size;
				int j = j_start + (temp % j_size);
				temp /= j_size;
				int k = k_start + temp;
				if (VElems(i, j, k) == porValue)
				{
					const int *pijk = &voxls(i, j, k);
					if (bgn <= voxls.v_i(1, pijk) && voxls.v_i(1, pijk) <= lst)
					{
						VElems(i, j, k) = voxls.v_i(1, pijk);
						++local_nChanges;
					}
					else if (bgn <= voxls.v_i(-1, pijk) && voxls.v_i(-1, pijk) <= lst)
					{
						VElems(i, j, k) = voxls.v_i(-1, pijk);
						++local_nChanges;
					}
					else if (bgn <= voxls.v_j(1, pijk) && voxls.v_j(1, pijk) <= lst)
					{
						VElems(i, j, k) = voxls.v_j(1, pijk);
						++local_nChanges;
					}
					else if (bgn <= voxls.v_j(-1, pijk) && voxls.v_j(-1, pijk) <= lst)
					{
						VElems(i, j, k) = voxls.v_j(-1, pijk);
						++local_nChanges;
					}
					else if (bgn <= voxls.v_k(1, pijk) && voxls.v_k(1, pijk) <= lst)
					{
						VElems(i, j, k) = voxls.v_k(1, pijk);
						++local_nChanges;
					}
					else if (bgn <= voxls.v_k(-1, pijk) && voxls.v_k(-1, pijk) <= lst)
					{
						VElems(i, j, k) = voxls.v_k(-1, pijk);
						++local_nChanges;
					}
				}
			}
			return local_nChanges;
		});
	for (auto &future : nChanges_futures)
	{
		nChanges += future.get();
	}
	cout << "  ngrowPors:" << nChanges << "  ";
}

void retreatPoresMedian(const inputDataNE &cg, voxelField<int> &VElems, long bgn, long lst,
						const vector<poreNE *> &poreIs, long unassigned)
{
	auto &pool = GlobalThreadPool::get();
	const voxelField<int> voxls = VElems;
	const int nz = VElems.nz(), ny = VElems.ny();
	// 内部区域
	const int j_start = 1, j_end = ny - 1;
	const int k_start = 1, k_end = nz - 1;
	const int j_size = j_end - j_start;
	const int k_size = k_end - k_start;
	const size_t total_kj = static_cast<size_t>(k_size * j_size); // collapse(2)
	auto futures = pool.submit_blocks(
		0, total_kj,
		[&](const size_t start_idx, const size_t end_idx) -> size_t
		{
			size_t local_nChanges = 0;

			for (size_t idx = start_idx; idx < end_idx; ++idx)
			{
				// 一维 idx -> 二维 (k, j)
				int j = j_start + (idx % j_size);
				int k = k_start + (idx / j_size);
				const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
				for (short ix = 0; ix < s.cnt; ++ix)
				{
					for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i)
					{
						const int *pijk = &voxls(i, j, k);
						long pID = *pijk;
						short nSameID = 0;
						short nDifferentID = 0;

						if (pID >= bgn && pID <= lst)
						{
							if (voxls.v_i(-1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_i(-1, pijk) && lst >= voxls.v_i(-1, pijk))
								nDifferentID++;
							if (voxls.v_i(1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk))
								nDifferentID++;
							if (voxls.v_j(-1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_j(-1, pijk) && lst >= voxls.v_j(-1, pijk))
								nDifferentID++;
							if (voxls.v_j(1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk))
								nDifferentID++;
							if (voxls.v_k(-1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_k(-1, pijk) && lst >= voxls.v_k(-1, pijk))
								nDifferentID++;
							if (voxls.v_k(1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk))
								nDifferentID++;

							if (nDifferentID > 0 && nSameID > 0)
							{
								VElems(i, j, k) = unassigned;
								++local_nChanges;
							}
						}
					}
				}
			}
			return local_nChanges;
		});
	size_t nChanges(0);
	for (auto &future : futures)
	{
		nChanges += future.get();
	}
	(cout << "  nRetreat:" << nChanges).flush();
}

void growPoresMedStrict(const inputDataNE &cg, voxelField<int> &VElems, long bgn, long lst,
						const vector<poreNE *> &poreIs, long rawValue)
{

	auto &pool = GlobalThreadPool::get();
	const voxelField<int> voxls = VElems;
	const int nz = VElems.nz(), ny = VElems.ny();
	const int j_start = 1, j_end = ny - 1;
	const int k_start = 1, k_end = nz - 1;
	const int j_size = j_end - j_start;
	const int k_size = k_end - k_start;
	const size_t total_iterations = static_cast<size_t>(k_size * j_size);
	auto nChanges_futures = pool.submit_blocks(
		0, total_iterations,
		[&](const size_t start_idx, const size_t end_idx) -> size_t
		{
			size_t local_nChanges = 0;
			for (size_t idx = start_idx; idx < end_idx; ++idx)
			{
				int j = j_start + (idx % j_size);
				int k = k_start + (idx / j_size);
				const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
				for (short ix = 0; ix < s.cnt; ++ix)
				{
					for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i)
					{
						long pID = voxls(i, j, k);
						const int *pijk = &voxls(i, j, k);

						if (pID == rawValue)
						{
							float R = cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 1)->R;
							short nDifferentID = 0;
							if (bgn <= voxls.v_i(-1, pijk) && lst >= voxls.v_i(-1, pijk) && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R >= R)
								nDifferentID++;
							if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk) && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R >= R)
								nDifferentID++;
							if (bgn <= voxls.v_j(-1, pijk) && lst >= voxls.v_j(-1, pijk) && cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R >= R)
								nDifferentID++;
							if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk) && cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R >= R)
								nDifferentID++;
							if (bgn <= voxls.v_k(-1, pijk) && lst >= voxls.v_k(-1, pijk) && cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
								nDifferentID++;
							if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk) && cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
								nDifferentID++;

							if (nDifferentID >= 3)
							{
								map<int, short> neis;

								long
									neI = voxls.v_i(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_i(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_j(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_j(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_k(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_k(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);

								map<int, short>::iterator neitr = max_element(neis.begin(), neis.end(), mapComparer<int>());
								if (neitr->second >= 3)
								{
									++local_nChanges;
									VElems(i, j, k) = neitr->first;
								}
							}
						}
					}
				}
			}
			return local_nChanges;
		});
	size_t nChanges(0);
	for (auto &future : nChanges_futures)
	{
		nChanges += future.get();
	}

	cout << "  ngMedStrict:" << nChanges << " ";
}

void growPoresMedian(const inputDataNE &cg, voxelField<int> &VElems, long bgn, long lst,
					 const vector<poreNE *> &poreIs, long rawValue)
{
	auto &pool = GlobalThreadPool::get();
	const voxelField<int> voxls = VElems;
	const int nz = VElems.nz(), ny = VElems.ny();
	const int j_start = 1, j_end = ny - 1;
	const int k_start = 1, k_end = nz - 1;
	const int j_size = j_end - j_start;
	const int k_size = k_end - k_start;
	const size_t total_iterations = static_cast<size_t>(k_size * j_size);
	auto nChanges_futures = pool.submit_blocks(
		0, total_iterations,
		[&](const size_t start_idx, const size_t end_idx) -> size_t
		{
			size_t local_nChanges = 0;
			for (size_t idx = start_idx; idx < end_idx; ++idx)
			{
				int j = j_start + (idx % j_size);
				int k = k_start + (idx / j_size);
				const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
				for (short ix = 0; ix < s.cnt; ++ix)
				{
					for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i)
					{
						// voxel* v=s.s[ix].v(i-1);
						// if (v)

						const int *pijk = &voxls(i, j, k);
						long pID = *pijk;

						if (pID == rawValue)
						{
							float R = cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 1)->R;

							short nDifferentID = 0;
							if (bgn <= voxls.v_i(-1, pijk) && lst >= voxls.v_i(-1, pijk) && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R > R)
								nDifferentID++;
							if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk) && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R > R)
								nDifferentID++;
							if (bgn <= voxls.v_j(-1, pijk) && lst >= voxls.v_j(-1, pijk) && cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R > R)
								nDifferentID++;
							if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk) && cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R > R)
								nDifferentID++;
							if (bgn <= voxls.v_k(-1, pijk) && lst >= voxls.v_k(-1, pijk) && cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R > R)
								nDifferentID++;
							if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk) && cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R > R)
								nDifferentID++;

							if (nDifferentID >= 2)
							{
								map<int, short> neis;

								long
									neI = voxls.v_i(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_i(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_j(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_j(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_k(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_k(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R > R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);

								map<int, short>::iterator neitr = max_element(neis.begin(), neis.end(), mapComparer<int>());
								if (neitr->second >= 2)
								{
									++local_nChanges;
									VElems(i, j, k) = neitr->first;
								}
							}
						}
					}
				}
			}
			return local_nChanges;
		});
	size_t nChanges(0);
	for (auto &future : nChanges_futures)
	{
		nChanges += future.get();
	}
	cout << "  ngMedian:" << nChanges << " ";
}

void growPoresMedEqs(const inputDataNE &cg, voxelField<int> &VElems, long bgn, long lst,
					 const vector<poreNE *> &poreIs, long rawValue)
{
	auto &pool = GlobalThreadPool::get();
	const voxelField<int> voxls = VElems;
	const int nz = VElems.nz(), ny = VElems.ny();
	const int j_start = 1, j_end = ny - 1;
	const int k_start = 1, k_end = nz - 1;
	const int j_size = j_end - j_start;
	const int k_size = k_end - k_start;
	const size_t total_iterations = static_cast<size_t>(k_size * j_size);
	auto nChanges_futures = pool.submit_blocks(
		0, total_iterations,
		[&](const size_t start_idx, const size_t end_idx) -> size_t
		{
			size_t local_nChanges = 0;
			for (size_t idx = start_idx; idx < end_idx; ++idx)
			{
				int j = j_start + (idx % j_size);
				int k = k_start + (idx / j_size);
				const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
				for (short ix = 0; ix < s.cnt; ++ix)
				{
					for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i)
					{
						// voxel* v=s.s[ix].v(i-1);
						// if (v)

						const int *pijk = &voxls(i, j, k);
						long pID = *pijk;

						if (pID == rawValue)
						{
							float R = cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 1)->R;

							short nDifferentID = 0;
							if (bgn <= voxls.v_i(-1, pijk) && lst >= voxls.v_i(-1, pijk) && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R >= R)
								nDifferentID++;
							if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk) && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R >= R)
								nDifferentID++;
							if (bgn <= voxls.v_j(-1, pijk) && lst >= voxls.v_j(-1, pijk) && cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R >= R)
								nDifferentID++;
							if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk) && cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R >= R)
								nDifferentID++;
							if (bgn <= voxls.v_k(-1, pijk) && lst >= voxls.v_k(-1, pijk) && cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
								nDifferentID++;
							if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk) && cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
								nDifferentID++;

							if (nDifferentID >= 2)
							{
								map<int, short> neis;

								long neI = voxls.v_i(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 2)->R >= R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_i(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i)->R >= R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_j(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + (j - 2)].vxl(i - 1)->R >= R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_j(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 1) * cg.ny + j].vxl(i - 1)->R >= R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_k(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[(k - 2) * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_k(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI && cg.segs_[k * cg.ny + (j - 1)].vxl(i - 1)->R >= R)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);

								map<int, short>::iterator neitr = max_element(neis.begin(), neis.end(), mapComparer<int>());
								if (neitr->second >= 2)
								{
									++local_nChanges;
									VElems(i, j, k) = neitr->first;
								}
							}
						}
					}
				}
			}
			return local_nChanges;
		});
	size_t nChanges(0);
	for (auto &future : nChanges_futures)
	{
		nChanges += future.get();
	}

	cout << "  ngMedEqs:" << nChanges << "  ";
}

void growPoresMedEqsLoose(const inputDataNE &cg, voxelField<int> &VElems, long bgn, long lst,
						  const vector<poreNE *> &poreIs, long rawValue)
{
	auto &pool = GlobalThreadPool::get();
	const voxelField<int> voxls = VElems;
	const int nz = VElems.nz(), ny = VElems.ny();
	const int j_start = 1, j_end = ny - 1;
	const int k_start = 1, k_end = nz - 1;
	const int j_size = j_end - j_start;
	const int k_size = k_end - k_start;
	const size_t total_iterations = static_cast<size_t>(k_size * j_size);
	auto nChanges_futures = pool.submit_blocks(
		0, total_iterations,
		[&](const size_t start_idx, const size_t end_idx) -> size_t
		{
			size_t local_nChanges = 0;
			for (size_t idx = start_idx; idx < end_idx; ++idx)
			{
				int j = j_start + (idx % j_size);
				int k = k_start + (idx / j_size);
				const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
				for (short ix = 0; ix < s.cnt; ++ix)
				{
					for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i)
					{
						// voxel* v=s.s[ix].v(i-1);
						// if (v)

						const int *pijk = &voxls(i, j, k);
						long pID = *pijk;

						if (pID == rawValue)
						{
							//			 float R = cg.segs_[(k - 1) * cg.ny + (j - 1)].vxl(i - 1)->R;

							short nDifferentID = 0;
							if (bgn <= voxls.v_i(-1, pijk) && lst >= voxls.v_i(-1, pijk))
								nDifferentID++;
							if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk))
								nDifferentID++;
							if (bgn <= voxls.v_j(-1, pijk) && lst >= voxls.v_j(-1, pijk))
								nDifferentID++;
							if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk))
								nDifferentID++;
							if (bgn <= voxls.v_k(-1, pijk) && lst >= voxls.v_k(-1, pijk))
								nDifferentID++;
							if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk))
								nDifferentID++;

							if (nDifferentID >= 2)
							{
								map<int, short> neis;

								long
									neI = voxls.v_i(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_i(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_j(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_j(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_k(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_k(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);

								map<int, short>::iterator neitr = max_element(neis.begin(), neis.end(), mapComparer<int>());
								if (neitr->second >= 2)
								{
									++local_nChanges;
									VElems(i, j, k) = neitr->first;
								}
							}
						}
					}
				}
			}
			return local_nChanges;
		});
	size_t nChanges(0);
	for (auto &future : nChanges_futures)
	{
		nChanges += future.get();
	}

	cout << "  ngMedLoose:" << nChanges << "  ";
}

void medianElem(const inputDataNE &cg, voxelField<int> &VElems, long bgn, long lst,
				const vector<poreNE *> &poreIs)
{
	auto &pool = GlobalThreadPool::get();
	const voxelField<int> voxls = VElems;
	const int nz = VElems.nz(), ny = VElems.ny();
	const int j_start = 1, j_end = ny - 1;
	const int k_start = 1, k_end = nz - 1;
	const int j_size = j_end - j_start;
	const int k_size = k_end - k_start;
	const size_t total_iterations = static_cast<size_t>(k_size * j_size);
	auto nChanges_futures = pool.submit_blocks(
		0, total_iterations,
		[&](const size_t start_idx, const size_t end_idx) -> size_t
		{
			size_t local_nChanges = 0;
			for (size_t idx = start_idx; idx < end_idx; ++idx)
			{
				int j = j_start + (idx % j_size);
				int k = k_start + (idx / j_size);
				const segments &s = cg.segs_[(k - 1) * cg.ny + (j - 1)];
				for (short ix = 0; ix < s.cnt; ++ix)
				{
					for (short i = s.s[ix].start + 1; i <= s.s[ix + 1].start; ++i)
					{
						const int *pijk = &voxls(i, j, k);
						long pID = *pijk;
						short nSameID = 0;
						short nDifferentID = 0;

						if (pID >= bgn && pID <= lst)
						{
							if (voxls.v_i(-1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_i(-1, pijk) && lst >= voxls.v_i(-1, pijk))
								nDifferentID++;
							if (voxls.v_i(1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_i(1, pijk) && lst >= voxls.v_i(1, pijk))
								nDifferentID++;
							if (voxls.v_j(-1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_j(-1, pijk) && lst >= voxls.v_j(-1, pijk))
								nDifferentID++;
							if (voxls.v_j(1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_j(1, pijk) && lst >= voxls.v_j(1, pijk))
								nDifferentID++;
							if (voxls.v_k(-1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_k(-1, pijk) && lst >= voxls.v_k(-1, pijk))
								nDifferentID++;
							if (voxls.v_k(1, pijk) == pID)
								nSameID++;
							else if (bgn <= voxls.v_k(1, pijk) && lst >= voxls.v_k(1, pijk))
								nDifferentID++;

							if (nDifferentID > nSameID)
							{
								map<int, short> neis;

								long
									neI = voxls.v_i(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_i(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_j(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_j(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_k(-1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								neI = voxls.v_k(1, pijk);
								if (neI != pID && bgn <= neI && lst >= neI)
									++(neis.insert(pair<int, short>(neI, 0)).first->second);
								for (map<int, short>::iterator neitr = neis.begin(); neitr != neis.end(); ++neitr)
								{
									if (neitr->second > nSameID)
									{
										++local_nChanges;
										VElems(i, j, k) = neitr->first;
										nSameID = neitr->second;
									}
								}
							}
						}
					}
				}
			}
			return local_nChanges;
		});
	size_t nChanges(0);
	for (auto &future : nChanges_futures)
	{
		nChanges += future.get();
	}

	cout << "  nMedian:" << nChanges << " ";
}

#endif
