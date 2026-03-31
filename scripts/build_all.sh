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

script_targets=( "${PROJ_DIR}/scripts/"*.cc )
datashaping_targets=( "${PROJ_DIR}/datashaping/"*.cpp )
check_targets=( "${PROJ_DIR}/check/eventdisplay.cc" )

targets=(
  "${script_targets[@]}"
  "${datashaping_targets[@]}"
  "${check_targets[@]}"
)

if (( ${#targets[@]} == 0 )); then
  echo "[error] no build targets found"
  exit 1
fi

# scripts/*.cc は src/*.cc とリンクしてライブラリ実装を使う
srcs=( "${PROJ_DIR}/src/"*.cc )

echo "[info] project : ${PROJ_DIR}"
echo "[info] build   : ${BUILD_DIR}"
echo "[info] compiler: ${CXX}"
echo "[info] scripts : ${#script_targets[@]}"
echo "[info] datashaping: ${#datashaping_targets[@]}"
echo "[info] check   : ${#check_targets[@]}"
echo "[info] targets : ${#targets[@]}"

n_ok=0
n_fail=0

for s in "${targets[@]}"; do
  base="$(basename "${s}")"
  name="${base%.*}"
  out="${BUILD_DIR}/${name}"
  extra_srcs=()

  case "${s}" in
    "${PROJ_DIR}/scripts/"*)
      extra_srcs=( "${srcs[@]}" )
      ;;
  esac

  echo "------------------------------------------------------------"
  echo "[build] ${base} -> build/${name}"

  # ここで 1つでも失敗したら止める（set -e）
  if (( ${#extra_srcs[@]} > 0 )); then
    "${CXX}" ${CXXFLAGS_BASE} ${INCLUDE_FLAGS} \
      ${ROOT_CFLAGS} \
      -o "${out}" \
      "${s}" \
      "${extra_srcs[@]}" \
      ${ROOT_LIBS}
  else
    "${CXX}" ${CXXFLAGS_BASE} ${INCLUDE_FLAGS} \
      ${ROOT_CFLAGS} \
      -o "${out}" \
      "${s}" \
      ${ROOT_LIBS}
  fi

  echo "[ok] build/${name}"
  n_ok=$((n_ok + 1))
done

echo "------------------------------------------------------------"
echo "[done] ok=${n_ok}  fail=${n_fail}"
