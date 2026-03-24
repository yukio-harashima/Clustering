# Finite Fault Clustering & Plotting Pipeline
PDTIの出力結果を用いて、階層型クラスタリングによるるいじの抽出と、GMT6による多様な可視化を行うための解析ツール群です。
- `1_clustering/`: PDTIの結果からスナップショットを生成し、クラスタリングを行うPythonスクリプトおよび高速化モジュール。
  
- `2_plotting/`: クラスタリング結果やインバージョン結果をGMTで描画・可視化するためのBashスクリプトおよびデータ整形用補助コード。
