# '''Map arg parser arguments to API arguments'''
# def map_args_filter_keep_reads(ap):
#     '''Filter keep reads mappings'''
#     opts = []
#     for field in [ap.r, ap.x, ap.rT, ap.xT, ap.lineagefile]:
#         try:
#             opts.append(field)
#         except AttributeError:
#             opts.append("")
#     return {
#         "input_file":  ap.i,
#         "kraken": ap.k,
#         "RetainIds": opts[0],
#         "ExcludeIds": opts[1],
#         "RetainNames": opts[2],
#         "ExcludeNames": opts[3],
#         "LineageFile": opts[4]
#     }

# def map_args_analysis(ap):
#     '''Analysis mappings'''
#     opts = []
#     for field in [ap.keepdups, ap.clin, ap.depth_inf]:
#         try:
#             opts.append(field)
#         except AttributeError:
#             opts.append("")
#     return {
#         "ExpName": ap.b,
#         "input_file":  ap.i,
#         "ExpDir": ap.o + "/",
#         "Probes": ap.p,
#         "Samples": ap.samples,
#         "KeepDups": opts[0],
#         "Clin": opts[1],
#         "DepthInf": opts[2]
#     }
