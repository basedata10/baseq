import base64, os, json

def build_result(process_path, outpath):
    from baseq.utils.buildresult import pack_web_datas
    print("[info] Build datas for visualization")
    data = pack_web_datas()
    data.add_image("genome_50k", os.path.join(process_path, "CNV_plot.50K.png"))
    data.add_image("genome_200k", os.path.join(process_path, "CNV_plot.200K.png"))
    data.add_image("genome_1M", os.path.join(process_path, "CNV_plot.1M.png"))
    data.add_image("GCBias", os.path.join(process_path, "GC_counts.png"))
    data.data["total"] = 5000000
    data.data["mapped"] = 4000000
    data.data["mapping_ratio"] = 0.8
    data.data["MAD"] = 0.5

    data.write_file(outpath)