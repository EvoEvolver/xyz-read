#!/bin/sh
set -eu

repo="EvoEvolver/xyz-read"
version="${XYZ_READ_VERSION:-latest}"
install_dir="${XYZ_READ_INSTALL_DIR:-${HOME}/.local/bin}"

case "$(uname -s)" in
    Linux) os="unknown-linux-musl" ;;
    Darwin) os="apple-darwin" ;;
    *)
        echo "xyz-read: unsupported operating system: $(uname -s)" >&2
        exit 1
        ;;
esac

case "$(uname -m)" in
    x86_64 | amd64) arch="x86_64" ;;
    arm64 | aarch64) arch="aarch64" ;;
    *)
        echo "xyz-read: unsupported CPU architecture: $(uname -m)" >&2
        exit 1
        ;;
esac

target="${arch}-${os}"
asset="xyz-read-${target}.tar.gz"
if [ "$version" = "latest" ]; then
    base_url="https://github.com/${repo}/releases/latest/download"
else
    case "$version" in
        v*) tag="$version" ;;
        *) tag="v${version}" ;;
    esac
    base_url="https://github.com/${repo}/releases/download/${tag}"
fi

temp_dir=$(mktemp -d "${TMPDIR:-/tmp}/xyz-read-install.XXXXXX")
cleanup() {
    rm -rf "$temp_dir"
}
trap cleanup EXIT HUP INT TERM

download() {
    url="$1"
    output="$2"
    if command -v curl >/dev/null 2>&1; then
        curl --fail --location --silent --show-error "$url" --output "$output"
    elif command -v wget >/dev/null 2>&1; then
        wget --quiet "$url" --output-document="$output"
    else
        echo "xyz-read: curl or wget is required" >&2
        exit 1
    fi
}

echo "Downloading xyz-read ${version} for ${target}..."
download "${base_url}/${asset}" "${temp_dir}/${asset}"
download "${base_url}/${asset}.sha256" "${temp_dir}/${asset}.sha256"

expected=$(awk 'NR == 1 { print $1 }' "${temp_dir}/${asset}.sha256")
if command -v sha256sum >/dev/null 2>&1; then
    actual=$(sha256sum "${temp_dir}/${asset}" | awk '{ print $1 }')
elif command -v shasum >/dev/null 2>&1; then
    actual=$(shasum -a 256 "${temp_dir}/${asset}" | awk '{ print $1 }')
else
    echo "xyz-read: sha256sum or shasum is required to verify the download" >&2
    exit 1
fi

if [ -z "$expected" ] || [ "$expected" != "$actual" ]; then
    echo "xyz-read: SHA-256 verification failed for ${asset}" >&2
    exit 1
fi

tar -xzf "${temp_dir}/${asset}" -C "$temp_dir"
if [ ! -f "${temp_dir}/xyz-read" ]; then
    echo "xyz-read: release archive did not contain the xyz-read binary" >&2
    exit 1
fi

mkdir -p "$install_dir"
if command -v install >/dev/null 2>&1; then
    install -m 0755 "${temp_dir}/xyz-read" "${install_dir}/xyz-read"
else
    cp "${temp_dir}/xyz-read" "${install_dir}/xyz-read"
    chmod 0755 "${install_dir}/xyz-read"
fi

echo "Installed xyz-read to ${install_dir}/xyz-read"
case ":${PATH}:" in
    *":${install_dir}:"*) ;;
    *) echo "Add ${install_dir} to PATH to run xyz-read from your shell." ;;
esac
