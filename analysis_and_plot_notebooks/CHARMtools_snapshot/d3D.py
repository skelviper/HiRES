
import cooler
import pandas as pd

def ratings_norm(df):
    min_value = df["distance"].min()
    max_value = df["distance"].max()
    MIN = 0
    MAX = 1
    df["distance_adj"] = df["distance"].apply(
        lambda x: MIN + (MAX - MIN) / (1 + max_value - min_value) * (x - min_value))

    return df

#def read_file(filepath, bintable, type="center", chrom="chr1"):
def read_file(filepath,bintable,type="center",chrom="chr18"): #,start=None,end=None):

    c = cooler.Cooler(filepath)  # 用cooler读入

    # 读入并fetch数据量减小
    mat = c.pixels(join=True).fetch(region=(chrom + 'mat'))
    pat = c.pixels(join=True).fetch(region=(chrom + 'pat'))
    # else:
    #     mat = c.pixels(join=True).fetch(region=(chrom + 'mat', start, end))
    #     pat = c.pixels(join=True).fetch(region=(chrom + 'pat', start, end))
    example = pd.concat([mat, pat], axis=0).rename(columns={"count": "distance", "chrom1": "chrom"})

    example['pos1'] = ((example['start1'] + example['end1']) / 2).astype(int)
    example['pos2'] = ((example['start2'] + example['end2']) / 2).astype(int)

    # 如果要scale的话
    if (type == 'center'):
        example['band'] = example['pos2'] - example['pos1']  # 添加新的一列
        # 分组后进行标准化
        example = example.groupby(["band"]).apply(ratings_norm)
        example = example.rename(columns={"distance": "raw"})
        example = example.rename(columns={"distance_adj": "distance"})  # 重命名这一列

        print("example shape", example.shape)
        print("bintable shape", bintable.shape)

        example = pd.merge(target_bintable, example, how='left', on=["chrom", "pos1", "pos2"])[
            ['chrom', 'pos1', 'pos2', 'distance', 'raw']]
        print("example shape", example.shape)

        example = example[['chrom', 'pos1', 'pos2', 'distance', 'raw']]

        #     下面一部分拼接到右面
        lens = int(example.shape[0] / 2)
        example = pd.concat([example[0:lens].reset_index(drop=True), example[lens:]['distance'].reset_index(drop=True),
                             example[lens:]['raw'].reset_index(drop=True)], axis=1)

        example.columns = ["chrom", "pos1", "pos2", "distance", "raw", "distance2", "raw2"]
        example = example[['distance', 'distance2', "raw", "raw2"]]  # 只取distance1和distance2

    else:
        # 和bintable进行merge
        example = pd.merge(bintable, example, how='left', on=["chrom","pos1", "pos2"])[
            ['chrom', 'pos1', 'pos2', 'distance']]

        example = example[['chrom', 'pos1', 'pos2', 'distance']]

        #     下面一部分拼接到右面
        lens = int(example.shape[0] / 2)

        example = pd.concat([example[0:lens].reset_index(drop=True), example[lens:]['distance'].reset_index(drop=True)],
                            axis=1)
        example.columns = ["chrom", "pos1", "pos2", "distance", "distance2"]
        example = example[['distance', 'distance2']]  # 只取distance1和distance2

    return example


def load_mat(filepaths, bintable, threads=40, type="center", chrom="chr1(mat)",start = None,end = None):
    # 多进程
    from multiprocessing import Process, Manager, Pool
    import pandas as pd

    pool = Pool(processes=threads)  # 创建线程池
    pool_data_list = []

    # 创建空dataframe
    data = pd.DataFrame()

    for file_name in filepaths:
        pool_data_list.append(pool.apply_async(read_file, (file_name, bintable, type, chrom,start,end)))

    pool.close()
    pool.join()

    # 从data_list中取出数据,这里有问题了
    for pool_data in pool_data_list:
        data = pd.concat([data, pool_data.get()], axis=1)  # 按列拼接

    return data


def fdr(p_vals):
    from scipy.stats import rankdata
    # 排序
    ranked_p_values = rankdata(p_vals)
    # 次数
    lens = p_vals.shape[0]
    # 计算
    fdr = p_vals * lens / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr


# 对每一列的处理
def d3Dtest(data, method='t', m=100):
    import numpy as np

    from collections import Counter
    from scipy.stats import chi2_contingency, ttest_ind

    x = data[0:m]# celltype1的数据
    y = data[m:]# celltype2的数据

    if (x.notnull().sum() == 0 or y.notnull().sum() == 0):
        return 1

    if method == 't':
        r, p = ttest_ind(x, y, nan_policy='omit')  # 必须要加这个参数

    else:
        count_x = Counter(x)  # 统计x中有多少0和1
        count_y = Counter(y)  # 统计y中有多少0和1

        obs = np.array([[count_x[0], count_x[1]], [count_y[0], count_y[1]]])

        if (method == 'chi-square'):
            stat, p, dof, expected = chi2_contingency(obs)
        else:
            oddsr, p = fisher_exact(obs)

    return p


def d3D(mat1, mat2, bintable,adj_method="BH", fdr_thres=0.05, test_method="chi-square", resolution=20000,
        center="yes",threads=10):

    from pandarallel import pandarallel
    import pandas as pd
    # pandarallel.initialize(progress_bar=True)
    pandarallel.initialize(nb_workers=threads)

    if (center == "yes"):

        mat1_distance = mat1[mat1.index.str.contains("distance")]  # 细胞类型1
        mat2_distance = mat2[mat2.index.str.contains("distance")]  # 细胞类型2

        # 2个细胞类型的拼起来
        mat = pd.concat([mat1_distance, mat2_distance], axis=0)
        # 得到pvalue的列
        print("mat")
        print(mat.head())
        # p = pd.Series(mat.parallel_apply(d3Dtest,args=('t',mat1_distance.shape[0])))#内存不足时会报错
        p = pd.Series(mat.apply(d3Dtest, args=('t', mat1_distance.shape[0])))
        print("pv",p.value_counts())

        # 添加pvalue列
        sig = pd.concat([target_bintable.iloc[0:int(target_bintable.shape[0] / 2), :], p], axis=1)
        sig.columns = ['chrom', 'pos1', 'chrom2', 'pos2', 'pv']

        sig['chrom'] = sig['chrom'].str.replace('mat', '')
        sig['chrom2'] = sig['chrom2'].str.replace('mat', '')

        type1_raw = mat1[mat1.index.str.contains("raw")].mean()
        type2_raw = mat2[mat2.index.str.contains("raw")].mean()

        sig = pd.concat([sig, type1_raw, type2_raw], axis=1)
        sig.columns = ['chrom', 'pos', 'chrom2', 'pos2', 'pv', 'mat1_mean', 'mat2_mean']

    else:
        mat = pd.concat([mat1, mat2], axis=0)
        p = pd.Series(mat.parallel_apply(d3Dtest, method='t'))
        sig = pd.concat([bintable, p], axis=1)
        sig.columns = ['chrom', 'pos1', 'chrom2', 'pos2', 'pv']

    # 对p进行调整
    if (adj_method == "FDR"):
        sig['FDR'] = fdr(sig['pv'])

    else:
        sig['FDR'] = sig['pv']

    # 算一下那个差
    sig['diff'] = mat1[mat1.index.str.contains("distance")].mean() - mat2[mat2.index.str.contains("distance")].mean()

    # 过滤fdr符合要求的行
    #     sig = sig[sig['FDR'] < fdr_thres]#先不过滤了

    # 添加start和end
    # sig['start1'] = (sig['pos1'] - resolution / 2).astype(int)
    # sig['end1'] = (sig['pos1'] + resolution / 2).astype(int)
    # sig['start2'] = (sig['pos2'] - resolution / 2).astype(int)
    # sig['end2'] = (sig['pos2'] + resolution / 2).astype(int)

    # sig = sig.drop(["pos", "pos2"], axis=1)

    # 根据阈值进行过滤
    # sig = sig[sig['FDR'] < fdr_thres]

    return (sig)

# todo:
# filter chromosome here is not a good idea since we always need to get a specific region in the downstream analysis
# a test run need 1min29s
def read_bin_file(binpath,chrom=None,start=None,end=None):
    import pandas as pd
    bintable = pd.read_csv(binpath, sep='\t')
    #去掉括号
    bintable['chrom'] = bintable['chrom'].str.replace('(', "").str.replace(')', "")
    bintable['chrom2'] = bintable['chrom2'].str.replace('(', "").str.replace(')', "")
    bintable.columns = ['chrom', 'pos1', 'chrom2', 'pos2']
    #调整pos的偏差
    bintable['pos1'] = (bintable['pos1'] - 10000).astype(int)
    bintable['pos2'] = (bintable['pos2'] - 10000).astype(int)
    #去掉xy染色体
    bintable = bintable.loc[bintable['chrom'].str.contains('chr[0-9]+')]  # 筛选

    if(chrom):
        mat = chrom + "mat"
        pat = chrom + "pat"
        if(start):
            start = start - 10000
            end = end -10000
            return (bintable.query("(chrom=='{}' | chrom=='{}') & pos1 >= @start & pos2 <=@end".format(mat, pat)))
        else:
            return (bintable.query("chrom =='{}' | chrom=='{}'".format(mat, pat)))
    else:
        return(bintable)


if __name__ == '__main__':
    import os
    import pandas as pd
    import cooler
    import warnings
    warnings.filterwarnings('ignore')

    filenames1 = []
    filenames2 = []

    metadata = pd.read_csv("/shareb/zliu/analysis/hires_mouse_dev/metadata_qcpass.tsv", sep="\t")[
        ["cellname", "celltype", "cellcycle_threshold", "rmsd_20k"]]
    # 365
    cellnames1 = metadata.query('celltype == "early neurons" & cellcycle_threshold == "G0" & rmsd_20k < 1.5')[
        "cellname"].tolist()
    # 269
    cellnames2 = metadata.query('celltype == "mix late mesenchyme" & cellcycle_threshold == "G0" & rmsd_20k < 1.5')[
        "cellname"].tolist()

    object_path = "/shareb/zliu/analysis/hires_mouse_dev/figure3_related/d3d_res_analysis/obs_exp/cools_distance"

    for var in os.listdir(object_path):
        if (var.split(".")[0] in cellnames1):
            filenames1.append(os.path.join(object_path, var))
        if (var.split(".")[0] in cellnames2):
            filenames2.append(os.path.join(object_path, var))

    print("reading binfile")
    binpath = "/shareb/zliu/analysis/hires_mouse_dev/figure3_related/d3d_res_analysis/bintable.20k.2m.tsv"
    target_bintable = read_bin_file(binpath, "chr1")

    print(target_bintable.columns)
    print(target_bintable.head())

    print("loading mat")
    mat1 = load_mat(filenames1[0:2], target_bintable, threads=5, type="center", chrom="chr1").T
    mat2 = load_mat(filenames2[0:2], target_bintable, threads=5, type="center", chrom="chr1").T

    print("d3D")
    sig = d3D(mat1,mat2,target_bintable,adj_method="FDR",fdr_thres=0.05,test_method="t",resolution=20000,center="yes")

    sig.to_csv("sig.csv")