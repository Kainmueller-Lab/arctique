from album.runner.api import setup

path_download_macos = "https://www.blender.org/download/release/Blender4.0/blender-4.0.2-macos-arm64.dmg/"
path_download_linux = "https://www.blender.org/download/release/Blender4.0/blender-4.0.2-linux-x64.tar.xz/"
path_download_windows = "https://www.blender.org/download/release/Blender4.0/blender-4.0.2-windows-x64.msi/"
download_name_linux = "blender-4.0.2-linux-x64.tar.xz"
download_name_windows = "blender-4.0.2-linux-x64.zip"
download_name_macos = "blender-4.0.2-linux-x64.dmg"
path_run_linux = "blender-4.0.2-linux-x64/blender"
path_run_windows = "blender-4.0.2-linux-x64\\blender.exe"
path_run_macos = "Contents/MacOS/Blender"
rel_python_path_linux = "blender-4.0.2-linux-x64/4.0/python/bin/python3.10"
rel_python_path_windows = "blender-4.0.2-linux-x64\\4.0\\python\\bin\\python.exe"
rel_python_path_macos = ""  # todo: updateme


def install():
    import sys
    operation_system = sys.platform
    if operation_system == "linux":
        __install_linux()
    elif operation_system == "darwin":
        __install_macos()
    else:
        __install_windows()
    __install_package()

def __install_package():
    import subprocess
    from album.runner.api import get_app_path

    python = str(get_app_path().joinpath(_get_blender_python()))

    args = [python, "install_rhender.py"]
    subprocess.run(args)


def __install_linux():
    from album.runner.api import get_cache_path, extract_tar, get_app_path
    download_target = get_cache_path().joinpath(download_name_linux)
    download(path_download_linux, download_target)
    extract_tar(download_target, get_app_path())


def __install_macos():
    import dmglib
    from distutils.dir_util import copy_tree
    import os
    from album.runner.api import get_cache_path, get_app_path
    get_app_path().mkdir(exist_ok=True, parents=True)
    download_target = get_cache_path().joinpath(download_name_macos)
    download(path_download_macos, download_target)
    dmg = dmglib.DiskImage(download_target)

    if dmg.has_license_agreement():
        print('Cannot attach disk image.')
        return

    for mount_point in dmg.attach():
        for entry in os.listdir(mount_point):
            print('{} -- {}'.format(mount_point, entry))
        copy_tree(mount_point + os.sep + 'Blender.app', str(get_app_path()))

    dmg.detach()


def __install_windows():
    from album.runner.api import get_cache_path, get_app_path
    download_target = get_cache_path().joinpath(download_name_windows)
    download(path_download_windows, download_target)
    extract_zip(download_target, get_app_path())


def extract_zip(in_zip, out_dir):
    from pathlib import Path
    import zipfile
    from album.runner.album_logging import get_active_logger
    out_path = Path(out_dir)

    if not out_path.exists():
        out_path.mkdir(parents=True)

    get_active_logger().info(f"Extracting {in_zip} to {out_dir}...")

    with zipfile.ZipFile(in_zip, 'r') as zip_ref:
        zip_ref.extractall(out_dir)


def _get_session():
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3 import Retry
    s = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.5)

    adapter = HTTPAdapter(max_retries=retry)

    s.mount("http://", adapter)
    s.mount("https://", adapter)

    return s


def _request_get(url):
    """Get a response from a request to a resource url."""
    from album.ci.utils.zenodo_api import ResponseStatus
    with _get_session() as s:
        r = s.get(url, allow_redirects=True, stream=True)

        if r.status_code != ResponseStatus.OK.value:
            raise ConnectionError("Could not connect to resource %s!" % url)

        return r


def download(str_input, output):
    with _get_session() as s:
        r = s.get(str_input, allow_redirects=True, stream=True)
        with open(output, "wb") as out:
            out.write(r.content)


def _get_blender_executable():
    import sys
    operation_system = sys.platform
    if operation_system == "linux":
        return path_run_linux
    elif operation_system == "darwin":
        return path_run_macos
    else:
        return path_run_windows

def _get_blender_python():
    import sys
    operation_system = sys.platform
    if operation_system == "linux":
        return rel_python_path_linux
    elif operation_system == "darwin":
        return rel_python_path_macos
    else:
        return rel_python_path_windows

def run_blender_script(script, *params):
    import subprocess
    from album.runner.api import get_app_path, get_args
    blender_path = str(get_app_path().joinpath(_get_blender_executable()))

    if get_args().headless:
        args = [blender_path, "-b", "-d", "-noaudio", "--debug-gpu-force-workarounds", "-P", script, "--"]
    else:
        args = [blender_path, "-d", "-noaudio", "--debug-gpu-force-workarounds", "-P", script, "--"]
    args.extend(params)
    subprocess.run(args)


def run():
    from pathlib import Path
    from album.runner.api import get_args, get_package_path

    project = Path(get_args().out)
    output_path = project.joinpath("export", "render")

    if output_path:
        Path(output_path).parent.mkdir(exist_ok=True, parents=True)

    run_blender_script(str(Path(get_package_path()).joinpath("main.py").absolute()))


setup(
    group="mdc-berlin",
    name="rHEnder",
    version="0.1.0",
    album_api_version="0.5.3",
    solution_creators=['Jan Philipp Albrecht', 'Deborah Schmidt'],
    cite=[{
        "text": "Blender Online Community: Blender - a 3D modelling and rendering package (2018). Stichting Blender Foundation, Amsterdam.",
        "url": "http://www.blender.org"
    }],
    title="rHEnder: Create H&E stained images in blender",
    description="",
    covers=[{
        "description": "",
        "source": "cover.png"
    }],
    run=run,
    install=install,
    args=[
        {
            "name": "out",
            "type": "directory",
            "description": "The output directory",
            "required": True
        },
        {
            "name": "headless",
            "default": False,
            "type": "boolean",
            "description": "Run Blender in the background or open the window (default)."
        }
    ],
    dependencies={'environment_file': """channels:
  - defaults
dependencies:
  - python=3.10
  - pip
  - requests
  - pip:
    - dmglib
"""}
)
