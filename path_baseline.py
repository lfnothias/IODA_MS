import numpy as np


def ReadFile(infile_name):
    data = my_data = np.genfromtxt(infile_name, delimiter=",", skip_header=1)
    return data


def DataFilter(data, intensity, intensity_ratio):
    data = data[data[:, 4] != 0]  # remove samples with intensity = 0
    data = data[data[:, 4] >= intensity]  # remove samples with intensity < given
    data = data[
        data[:, 4] / (data[:, 3] + 1e-4) > intensity_ratio
    ]  # remove samples with intensity ratio < given
    return data


def PathGen(data, window_len, num_path, iso):
    start = min(data[:, 1])
    end = max(data[:, 1])
    path = []
    while start < end:
        curr_end = start + float(window_len)
        tmp_data = data[(data[:, 1] >= start) & (data[:, 1] < curr_end)]
        if len(tmp_data) != 0:
            ind = np.argsort(tmp_data[:, 4])
            ind = ind[::-1]
            tmp = []
            for i in range(num_path):
                if i >= len(ind):
                    break
                tmp.append(
                    (
                        tmp_data[ind[i], 0],
                        iso,
                        window_len,
                        start,
                        curr_end,
                        tmp_data[ind[i], 4],
                        tmp_data[ind[i], 1],
                        tmp_data[ind[i], 2],
                    )
                )
            path.append(tmp)
        start = curr_end
    return path


def WriteFile(outfile_name, path):
    text_file = open(outfile_name, "wt")
    for i in range(len(path[0])):
        n = text_file.write("path" + str(i) + "\t")
        for j in range(len(path)):
            if i > len(path[j]) - 1:
                continue
            n = text_file.write(
                "{:.4f}".format(path[j][i][0])
                + " "
                + "{:.4f}".format(path[j][i][1])
                + " "
                + "{:.4f}".format(path[j][i][2])
                + " "
                + "{:.4f}".format(path[j][i][3])
                + " "
                + "{:.4f}".format(path[j][i][4])
                + " "
                + "{:.4f}".format(path[j][i][5])
                + " "
                + "{:.4f}".format(path[j][i][6])
                + " "
                + str(path[j][i][7])
                + "\t"
            )
        n = text_file.write("\n")
    text_file.close()
