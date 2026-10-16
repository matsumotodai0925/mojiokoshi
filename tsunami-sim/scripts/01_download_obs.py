"""
検証用データの取得（2025-07-29 M8.8 カムチャッカ地震）。

  data/usgs_ffm/  : USGS 震源情報(us6000qw60)と有限断層モデル(complete_inversion.fsp ほか)
  data/dart/      : NOAA NDBC の DART 2025 年データ（15 分値＋イベント時 1 分値）
  data/tide_ioc/  : IOC Sea Level Station Monitoring Facility の検潮記録（1 分値）
                    ※日本の検潮所の値は気象庁から IOC に提供されているもの
使い方: python3 01_download_obs.py
"""
import os
import json
import time
import urllib.request
import common as C

EVENT = 'us6000qw60'
DART = [21416, 21419, 21413, 21418, 21414, 21415, 46408, 52402, 51407]
# 潮汐調和解析のため地震の前後数日をとる
IOC_START, IOC_END = '2025-07-26T00:00', '2025-08-02T00:00'


def get(url, fn, retries=4):
    for i in range(retries):
        try:
            urllib.request.urlretrieve(url, fn)
            return True
        except Exception as e:  # noqa: BLE001
            print('retry', url, e)
            time.sleep(2 ** (i + 1))
    return False


def usgs():
    d = os.path.join(C.DATA, 'usgs_ffm'); os.makedirs(d, exist_ok=True)
    ev = os.path.join(d, f'event_{EVENT}.json')
    get(f'https://earthquake.usgs.gov/fdsnws/event/1/query?eventid={EVENT}&format=geojson', ev)
    j = json.load(open(ev))
    ff = j['properties']['products']['finite-fault'][0]['contents']
    for k in ('complete_inversion.fsp', 'basic_inversion.param', 'FFM.geojson'):
        get(ff[k]['url'], os.path.join(d, k))


def dart():
    d = os.path.join(C.DATA, 'dart'); os.makedirs(d, exist_ok=True)
    for s in DART:
        get(f'https://www.ndbc.noaa.gov/data/historical/dart/{s}t2025.txt.gz', os.path.join(d, f'{s}t2025.txt.gz'))


def ioc():
    d = os.path.join(C.DATA, 'tide_ioc'); os.makedirs(d, exist_ok=True)
    for code in sorted(set(C.IOC_CODE.values())):
        url = ('https://www.ioc-sealevelmonitoring.org/service.php?query=data'
               f'&code={code}&timestart={IOC_START}&timestop={IOC_END}&format=json')
        get(url, os.path.join(d, f'{code}.json'))


if __name__ == '__main__':
    usgs(); dart(); ioc()
