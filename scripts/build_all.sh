#!/usr/bin/env bash
set -euo pipefail

# 一括ビルド用シェルスクリプト

# 全部ビルドする ./scripts/build_all.sh
# クリーン ./scripts/build_all.sh clean

cmd="${1:-build}"

ROOT_CFLAGS="$(root-config --cflags)"
ROOT_LIBS="$(root-config --libs)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJ_DIR}/build"

CXX="${CXX:-g++}"
CXXFLAGS_BASE="-O2 -std=c++17 -Wall -Wextra -pedantic"
INCLUDE_FLAGS="-I${PROJ_DIR}/include"

mkdir -p "${BUILD_DIR}"

if [[ "${cmd}" == "clean" ]]; then
  echo "[clean] remove build/*"
  rm -f "${BUILD_DIR}/"*
  exit 0
fi

shopt -s nullglob

scripts=( "${PROJ_DIR}/scripts/"*.cc )
if (( ${#scripts[@]} == 0 )); then
  echo "[error] no scripts found: ${PROJ_DIR}/scripts/*.cc"
  exit 1
fi

# src/*.cc をまとめてリンク（現状運用に合わせる）
srcs=( "${PROJ_DIR}/src/"*.cc )

echo "[info] project : ${PROJ_DIR}"
echo "[info] build   : ${BUILD_DIR}"
echo "[info] compiler: ${CXX}"
echo "[info] scripts : ${#scripts[@]}"

n_ok=0
n_fail=0

for s in "${scripts[@]}"; do
  base="$(basename "${s}")"
  name="${base%.cc}"
  out="${BUILD_DIR}/${name}"

  echo "------------------------------------------------------------"
  echo "[build] ${base} -> build/${name}"

  # ここで 1つでも失敗したら止める（set -e）
  "${CXX}" ${CXXFLAGS_BASE} ${INCLUDE_FLAGS} \
    ${ROOT_CFLAGS} \
    -o "${out}" \
    "${s}" \
    "${srcs[@]}" \
    ${ROOT_LIBS}

  echo "[ok] build/${name}"
  n_ok=$((n_ok + 1))
done

echo "------------------------------------------------------------"
echo "[done] ok=${n_ok}  fail=${n_fail}"
