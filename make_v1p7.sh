#!/usr/bin/env bash
# ----------------------------------------------------------------------
# File     : make_v1p7.sh  (run once from the bearing repo root)
# Version  : 1.0.0
# Date     : 2026-10-20
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# Builds the V1P 7-track inputs (8-track panel minus the 4C Trbv13-2 channel;
# decision Craig + Nina 2026-10-20):
#   workflow/config/v1_7track_panel.yaml
#   workflow/config/samples_v1p_7track.tsv   (last bw of each row dropped)
#   workflow/config/config_v1p7_prior1_l1.yaml (diff_stat l1, new outdir)
# Refuses to overwrite existing files. ASCII-only.
set -euo pipefail
C=workflow/config
for f in $C/v1_7track_panel.yaml $C/samples_v1p_7track.tsv $C/config_v1p7_prior1_l1.yaml; do
  [ -e "$f" ] && { echo "exists, not overwriting: $f" >&2; exit 1; }
done
[ -f $C/samples_v1p_8track.tsv ] || { echo "missing $C/samples_v1p_8track.tsv" >&2; exit 1; }

# 1. panel: v1_8track_panel.yaml without the last category (4C Trbv13-2)
cat > $C/v1_7track_panel.yaml <<'YAML'
# ----------------------------------------------------------------------
# File     : workflow/config/v1_7track_panel.yaml
# Version  : 1.0.0
# Date     : 2026-10-20
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# V1P-vs-DN seven-track panel = v1_8track_panel minus the 4C Trbv13-2 channel.
# The Trbv13-2 O/E 4C channel carried 4,196 of 4,853 V1P BH bins across the
# whole 4C reach (chr6 35-48 Mb) and was V1P-higher at the RC; treated as
# non-specific (DN reference maps merge conventional + Arima libraries, V1P is
# Arima only). Decision Craig + Nina 2026-10-20: report V1P on the RC viewpoint
# plus the linear tracks.
# Order is FIXED and must match the bw column in samples_v1p_7track.tsv:
#   RNAseq +, RNAseq -, CTCF, Cohesin, NIPBL, H3K27ac, 4C RC
# ASCII-only.
name: "V1 7-track panel"
genome: mm10
categories:
  - name: "RNAseq +"
    color: "#6495ed"
  - name: "RNAseq -"
    color: "#1a3a8f"
    negative_strand: true
  - name: CTCF
    color: "#ff2200"
  - name: Cohesin
    color: "#8b0000"
  - name: NIPBL
    color: "#00bcd4"
  - name: H3K27ac
    color: "#00e676"
  - name: "4C RC"
    color: "#ffab00"
YAML

# 2. sample sheet: drop the 8th (last) bw path in every data row; keep comments
awk -F'\t' 'BEGIN{OFS="\t"}
  /^#/ {print; next}
  $1=="sample" {for(i=1;i<=NF;i++) if($i=="bw") b=i; print; next}
  { n=split($b, a, ","); if(n!=8){print "row "$1": expected 8 bw, found "n > "/dev/stderr"; bad=1}
    s=a[1]; for(i=2;i<n;i++) s=s","a[i]; $b=s; print }
  END{exit bad}' $C/samples_v1p_8track.tsv > $C/samples_v1p_7track.tsv
sed -i '1i # 7-track sheet made from samples_v1p_8track.tsv by make_v1p7.sh 1.0.0 (2026-10-20): 4C Trbv13-2 bw dropped.' $C/samples_v1p_7track.tsv

# 3. config: copy of config_v1p_prior1.yaml with new sheet, panel, outdir, diff_stat l1
sed -e 's#^samples_sheet:.*#samples_sheet: "config/samples_v1p_7track.tsv"#' \
    -e 's#^categories:.*#categories:    "config/v1_7track_panel.yaml"#' \
    -e 's#^outdir:.*#outdir: "results_v1p7_prior1_l1"#' \
    $C/config_v1p_prior1.yaml > $C/config_v1p7_prior1_l1.yaml
grep -q '^diff_stat:' $C/config_v1p7_prior1_l1.yaml \
  && sed -i 's#^diff_stat:.*#diff_stat: l1#' $C/config_v1p7_prior1_l1.yaml \
  || printf '\n# V1P 7-track run (2026-10-20, Nina + Claude): non-cancelling statistic\ndiff_stat: l1\n' >> $C/config_v1p7_prior1_l1.yaml
sed -i '1i # config_v1p7_prior1_l1.yaml: generated 2026-10-20 by make_v1p7.sh 1.0.0 (K. E. Hayer + Claude) from config_v1p_prior1.yaml; 7-track panel, diff_stat l1.' $C/config_v1p7_prior1_l1.yaml

echo "== check"
grep -nE '^(samples_sheet|categories|outdir|diff_stat|prior_strength|n_perms|min_signal|pvalue_min_signal):' $C/config_v1p7_prior1_l1.yaml
awk -F'\t' '/^#/{next} $1=="sample"{for(i=1;i<=NF;i++) if($i=="bw") b=i; next} {n=split($b,a,","); print $1, n" bw, last: "a[n]}' $C/samples_v1p_7track.tsv
