# submission-requirements — Physical Review Letters × 成果物 1対1照合（2026-10-29）

| # | 要件 | 状態 | 成果物 / 備考 |
|---|---|---|---|
| 1 | REVTeX 4.2 ソース（`prl` クラス、reprint） | ✓ | `main.tex` |
| 2 | Letter 長 ≤ 3750 語相当 | ✓ | 3715（本文 2875 語 + 式 5 行 + 図 4 枚、`prl_check.py` 換算） |
| 3 | アブストラクト（1 主張、≤ 196 語） | ✓ | 156 語、"Here we show"、平文の数値 1 つ、含意で閉じる |
| 4 | 図（本文が参照する全図） | ✓ | `figures/fig1–4.pdf`（ベクトル、`make_figures.py` で JSON から再生成） |
| 5 | 参照文献 | ✓ | `refs.bib` 47 件、引用 42 件、2024 年以降 31%、未定義引用 0 |
| 6 | Supplemental Material | ✓ | `sm.tex` / `sm.pdf`（S1–S7、表 4） |
| 7 | カバーレター（推奨査読者、投稿歴） | ✓ | `cover_letter.tex/pdf`（査読者 5 名、本稿・関連原稿の未投稿を明記） |
| 8 | Data Availability 節 | ✓ | `main.tex` §Data availability（Zenodo DOI は投稿時に挿入） |
| 9 | Competing interests / funding / contributions | ✓ | `DECLARATIONS.md` |
| 10 | 投稿衛生（改訂経緯の痕跡なし） | ✓ | `grep -inE "corrected|earlier (revision|draft)|previously reported|referee|phase-1" main.tex` ヒット 0 |
| 11 | コンパイル健全性 | ✓ | バンドル内で main / sm / cover_letter / arxiv すべてエラー 0 |
| 12 | arXiv 版（単一ファイル、bbl 同梱） | ✓ | `arxiv/main_arxiv.tex`, `main_arxiv.bbl`, `SUBMISSION_NOTES.md` |
| 13 | Zenodo 寄託物 | ✓（未寄託） | `zenodo/`（README, CITATION.cff, submission.json, tarball）。寄託・公開は著者承認後 |

## 投稿前に著者が埋める項目
- [ ] arXiv 番号（`cover_letter.tex`, `arxiv/SUBMISSION_NOTES.md`）
- [x] Zenodo: 本番 draft 作成済み（deposition 23801927、concept DOI 10.5281/zenodo.23801926）。DOI は転記済み。**公開（publish）は未実行**: https://zenodo.org/deposit/23801927 で内容確認後に Publish
- [ ] Author contributions の最終文言（`DECLARATIONS.md`）
- [ ] 推奨査読者の所属メール（投稿フォームで要求される場合）
