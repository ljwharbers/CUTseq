from collections import Counter

out_file = open("E:\SciLife Lab\CUTseq\Cut_site_analysis/distanceCounts.txt", "wt")

with open("E:\SciLife Lab\CUTseq\Cut_site_analysis/rb_batch_Subset_cutsites_distributionDistancesTabs.txt", "rt") as in_file:
    for line in in_file:
        line = line.split("\t")
        cnt = Counter()
        for i in range(2, len(line)):
            cnt.update([line[i]])
        for key in cnt:
            out_file.write(str(line[1]) + "\t" + str(key) + "\t" + str(cnt[key]) + "\n")
            #print(line[1], key, cnt[key], "\n")