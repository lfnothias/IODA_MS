import math
import os
import sys

import pandas as pd


def split_features(input_filename: str, output_filename:str, sample_name:str, rt_window: float, num_splits: int):
    features = pd.read_csv(input_filename).fillna(0).sort_values('retention_time')
    splits = list([] for _ in range(num_splits))
    rt = 0
    while rt < features['retention_time'].max():
        features_rt_window = features[
            (features['retention_time'] > rt) &
            (features['retention_time'] <= rt + rt_window)].sort_values(
                sample_name, ascending=False)
        if len(features_rt_window) > 0:
            split_len = math.ceil(len(features_rt_window) / num_splits)
            for split in range(min(num_splits, split_len + 1)):
                start_i, stop_i = split * split_len, (split + 1) * split_len
                splits[split].append(features_rt_window.iloc[start_i:stop_i])
        rt += rt_window
    _, ext = os.path.splitext(input_filename)
    for i, split in enumerate(splits):
        pd.concat(split).to_csv(output_filename.replace(ext, f'_{i+1}{ext}'),
                                index=False)
    output_filenames = output_filename.replace(ext, f'_{i}{ext}')

if __name__ == '__main__':
    split_features(sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]))
