# グラフ出力規則

- 出力ファイルは、種別ごとに個別フォルダを作成した上でそこに保存すること
- barplotにlegendを埋め込まず、barplotとlegendを別ファイルとして出力すること。
- 節足動物種ごとに個別フォルダを作り、対応するbarplotとlegendを同じ種別フォルダへ出力すること。
- barplotとlegendはSVG・PNGの両形式で出力すること。
- scatter plotにもlegendを埋め込まず、scatter plotとlegendを別ファイルとして同じ種別フォルダへ出力すること。
- scatter plotの縦寸法はbarplotへ合わせず、scatterの可読性に応じて独立に設定すること。
- 特定の菌だけを強調するなどの個別対応では、`gen_barplot`、`gen_scatter`などの大元のファイルを編集しないこと。大元の処理を参考に、一時スクリプトをBashから実行して別名の派生ファイルを出力すること。
- 派生ファイルでは、依頼で明示されていない寸法やレイアウトを変更しないこと。
