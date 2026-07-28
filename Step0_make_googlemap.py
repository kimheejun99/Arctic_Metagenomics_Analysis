#!/usr/bin/env python3
"""
Step0_make_googlemap.py

현재 위치의 샘플명/위도/경도 txt 파일을 구글맵 HTML 시각화 파일로 변환한다.

입력 파일 형식 (tab 또는 공백 구분, 헤더 필요):
    Run	latitude	longitude
    SRR1234	71.29	-156.5
    ...

사용법:
    python3 Step0_make_googlemap.py sample_coords.txt --outdir . --api_key YOUR_KEY
"""

import argparse
import json
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("coords_file", help="샘플명/위도/경도 텍스트 파일")
    p.add_argument("--sample_col", default="Run")
    p.add_argument("--lat_col", default="latitude")
    p.add_argument("--lon_col", default="longitude")
    p.add_argument("--api_key", default="YOUR_GOOGLE_MAPS_API_KEY", help="구글맵 API 키")
    p.add_argument("--outdir", default=".")
    p.add_argument("--out_name", default="sample_map.html")
    return p.parse_args()


HTML_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
  <title>Sample Map</title>
  <meta charset="utf-8">
  <style>
    #map {{ height: 100vh; width: 100%; }}
  </style>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/OverlappingMarkerSpiderfier/1.0.3/oms.min.js"></script>
</head>
<body>
  <div id="map"></div>
  <script>
    const data = {data_json};

    function initMap() {{
      const center = {{ lat: {center_lat}, lng: {center_lng} }};
      const map = new google.maps.Map(document.getElementById("map"), {{ zoom: 3, center: center }});
      const infoWindow = new google.maps.InfoWindow();
      const oms = new OverlappingMarkerSpiderfier(map, {{ markersWontMove: true, markersWontHide: true, keepSpiderfied: true }});

      data.forEach(loc => {{
        const marker = new google.maps.Marker({{
          position: {{ lat: loc.lat, lng: loc.lng }},
          title: loc.id
        }});
        marker.addListener('click', () => {{
          infoWindow.setContent(`<strong>${{loc.id}}</strong>`);
          infoWindow.open(map, marker);
        }});
        oms.addMarker(marker);
      }});
      oms.addListener('click', function(marker) {{
        infoWindow.open(map, marker);
      }});
    }}
  </script>
  <script async defer src="https://maps.googleapis.com/maps/api/js?key={api_key}&callback=initMap"></script>
</body>
</html>
"""


def main():
    args = parse_args()
    df = pd.read_csv(args.coords_file, sep=None, engine="python")

    data = [{"id": str(r[args.sample_col]), "lat": float(r[args.lat_col]), "lng": float(r[args.lon_col])}
            for _, r in df.iterrows()]

    center_lat = sum(d["lat"] for d in data) / len(data)
    center_lng = sum(d["lng"] for d in data) / len(data)

    html = HTML_TEMPLATE.format(
        data_json=json.dumps(data, ensure_ascii=False),
        center_lat=center_lat, center_lng=center_lng, api_key=args.api_key
    )

    out_path = f"{args.outdir}/{args.out_name}"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[Step0] 완료: {out_path} ({len(data)}개 샘플)")


if __name__ == "__main__":
    main()
