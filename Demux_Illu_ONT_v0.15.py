#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import logging
import multiprocessing as mp
import os
import sys
from datetime import datetime
from collections import defaultdict

import pandas as pd
from Bio import SeqIO


# ---------------- Логирование ---------------- #
def setup_logging(log_file):
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


# ---------------- Загрузка баз ---------------- #
def load_bases(base_file, index_file):
    base_df = pd.read_excel(base_file)
    index_df = pd.read_excel(index_file)

    s5_map = {row["Forward"]: row["Name"] for _, row in base_df.iterrows() if str(row["Name"]).startswith("S5")}
    n7_map = {row["Forward"]: row["Name"] for _, row in base_df.iterrows() if str(row["Name"]).startswith("N7")}
    s5_map_r = {row["Reverse"]: row["Name"] for _, row in base_df.iterrows() if str(row["Name"]).startswith("S5")}
    n7_map_r = {row["Reverse"]: row["Name"] for _, row in base_df.iterrows() if str(row["Name"]).startswith("N7")}

    index_map = {}
    for _, row in index_df.iterrows():
        s5 = str(row["S5"]).strip()
        n7 = str(row["N7"]).strip()
        index_map[(s5, n7)] = str(row["Sample"]).strip()

    return s5_map, n7_map, s5_map_r, n7_map_r, index_map


# ---------------- Анализ одного рида ---------------- #
def analyze_read(record, s5_map, n7_map, s5_map_r, n7_map_r, index_map):
    seq = str(record.seq)
    header = record.id

    start40 = seq[:40]
    end40 = seq[-40:]

    scheme, s5, n7, sample = "No chain direction", None, None, None

    if "GATCTACAC" in start40:  # Forward
        scheme = "Forward"
        for s5_seq, s5_name in s5_map.items():
            if ("GATCTACAC" + s5_seq) in start40:
                s5 = s5_name
                break
        for n7_seq, n7_name in n7_map.items():
            if (n7_seq + "ATCTCGTATGCC") in end40:
                n7 = n7_name
                break
    elif "CATACGAGAT" in start40:  # Reverse
        scheme = "Reverse"
        for n7_seq, n7_name in n7_map_r.items():
            if ("CATACGAGAT" + n7_seq) in start40:
                n7 = n7_name
                break
        for s5_seq, s5_name in s5_map_r.items():
            if (s5_seq + "GTGTAGATC") in end40:
                s5 = s5_name
                break

    if s5 and n7 and (s5, n7) in index_map:
        sample = index_map[(s5, n7)]

    return {
        "ReadID": header,
        "Scheme": scheme,
        "S5": s5,
        "N7": n7,
        "Sample": sample,
        "Record": record,
    }


# ---------------- Обработка результатов батча ---------------- #
def handle_results(results, sample_handles, unmultiplexed_handle, scheme_csv,
                   total, identified, unmapped, sample_stats):
    with open(scheme_csv, "a") as f:
        for r in results:
            total += 1
            if r["Sample"]:
                identified += 1
                if r["Sample"] not in sample_handles:
                    out_path = os.path.join(os.path.dirname(scheme_csv), f"{r['Sample']}.fastq.gz")
                    sample_handles[r["Sample"]] = gzip.open(out_path, "wt")
                SeqIO.write(r["Record"], sample_handles[r["Sample"]], "fastq")
                sample_stats[(r["Sample"], r["S5"], r["N7"])] += 1
            else:
                unmapped += 1
                SeqIO.write(r["Record"], unmultiplexed_handle, "fastq")
            f.write(f"{r['ReadID']},{r['Scheme']},{r['S5']},{r['N7']},{r['Sample']}\n")
    return total, identified, unmapped


# ---------------- Основная обработка fastq ---------------- #
def process_fastq(input_file, base_file, index_file, output_dir, threads, batch_size, max_excel_rows):
    s5_map, n7_map, s5_map_r, n7_map_r, index_map = load_bases(base_file, index_file)

    logging.info("Загрузка данных завершена")
    logging.info(f"Всего индексов S5: {len(s5_map) + len(s5_map_r)}, N7: {len(n7_map) + len(n7_map_r)}")

    if input_file.endswith(".gz"):
        handle = gzip.open(input_file, "rt")
    else:
        handle = open(input_file, "r")

    scheme_csv = os.path.join(output_dir, "reads_scheme.csv")
    with open(scheme_csv, "w") as f:
        f.write("ReadID,Scheme,S5,N7,Sample\n")

    sample_handles = {}
    unmultiplexed_path = os.path.join(output_dir, "unmultiplexed.fastq.gz")
    unmultiplexed_handle = gzip.open(unmultiplexed_path, "wt")

    total = 0
    identified = 0
    unmapped = 0
    sample_stats = defaultdict(int)

    batch = []
    pool = mp.Pool(threads)

    for record in SeqIO.parse(handle, "fastq"):
        batch.append(record)
        if len(batch) >= batch_size:
            results = pool.starmap(
                analyze_read,
                [(r, s5_map, n7_map, s5_map_r, n7_map_r, index_map) for r in batch],
            )
            total, identified, unmapped = handle_results(
                results, sample_handles, unmultiplexed_handle, scheme_csv,
                total, identified, unmapped, sample_stats
            )
            batch = []

    if batch:
        results = pool.starmap(
            analyze_read,
            [(r, s5_map, n7_map, s5_map_r, n7_map_r, index_map) for r in batch],
        )
        total, identified, unmapped = handle_results(
            results, sample_handles, unmultiplexed_handle, scheme_csv,
            total, identified, unmapped, sample_stats
        )

    pool.close()
    pool.join()
    handle.close()

    for h in sample_handles.values():
        h.close()
    unmultiplexed_handle.close()

    # ---------------- Сохранение схемы ---------------- #
    scheme_xlsx = os.path.join(output_dir, "reads_scheme.xlsx")
    df = pd.read_csv(scheme_csv)
    if len(df) > max_excel_rows:
        logging.warning(f"Количество строк ({len(df)}) превышает max_excel_rows={max_excel_rows}. Сохраняем только первые {max_excel_rows} строк в Excel")
        df.head(max_excel_rows).to_excel(scheme_xlsx, index=False)
    else:
        df.to_excel(scheme_xlsx, index=False)

    # ---------------- Сводка по образцам ---------------- #
    summary_csv = os.path.join(output_dir, "sample_summary.csv")
    summary_xlsx = os.path.join(output_dir, "sample_summary.xlsx")
    summary_data = [
        {"Sample": s, "S5": s5, "N7": n7, "Reads": c}
        for (s, s5, n7), c in sample_stats.items()
    ]
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(summary_csv, index=False)
    summary_df.to_excel(summary_xlsx, index=False)

    return total, identified, unmapped, scheme_csv, scheme_xlsx, summary_csv, summary_xlsx


# ---------------- main ---------------- #
def main():
    parser = argparse.ArgumentParser(description="Demultiplex FASTQ по индексам (streaming, low memory)")
    parser.add_argument("-i", "--input", required=True, help="Входной fastq(.gz)")
    parser.add_argument("-b", "--base", required=True, help="Файл base.xlsx")
    parser.add_argument("-x", "--index", required=True, help="Файл Index.xlsx")
    parser.add_argument("-o", "--output", required=True, help="Папка для результатов")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Количество потоков")
    parser.add_argument("--batch", type=int, default=100000, help="Размер батча (по умолчанию 100000)")
    parser.add_argument("--max-excel-rows", type=int, default=1000000, help="Максимальное число строк для Excel")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    log_file = os.path.join(args.output, "demultiplex.log")
    setup_logging(log_file)

    logging.info("Запуск скрипта демультиплексирования (версия 0.15, streaming)")
    logging.info(f"Команда: {' '.join(sys.argv)}")

    start = datetime.now()

    total, identified, unmapped, scheme_csv, scheme_xlsx, summary_csv, summary_xlsx = process_fastq(
        args.input, args.base, args.index, args.output, args.threads, args.batch, args.max_excel_rows
    )

    end = datetime.now()
    elapsed = (end - start).total_seconds()

    logging.info(f"Всего ридов: {total}")
    logging.info(f"Идентифицировано: {identified}")
    logging.info(f"Неидентифицировано: {unmapped}")
    logging.info(f"Файл схемы CSV: {scheme_csv}")
    logging.info(f"Файл схемы Excel: {scheme_xlsx}")
    logging.info(f"Сводка по образцам CSV: {summary_csv}")
    logging.info(f"Сводка по образцам Excel: {summary_xlsx}")
    logging.info(f"Время работы: {elapsed:.2f} секунд")


if __name__ == "__main__":
    main()
