#!/bin/bash
#$ -S /bin/bash

# This script is used to perform cleanup of the ASHS output by trimming slices that have
# only a few ERC/PRC voxels (less than 25% of the median number). It also generates some
# screenshots that are now largely integrated into the fast-ashs branch of ASHS and so
# are largely redundant

# The script assumes that there is a subdirectory called ashs that contains directories 
# corresponding to a set of subjects. Each subject's directory coincides with the subject's
# id. There is also a list of subject ids (separated by space) stored in cleanup_subj.txt

# (!) Path to the ImageMagick library
IMROOT=~/bin/imagemagick
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$IMROOT/lib
export PATH=$PATH:$IMROOT/bin

mkdir -p cleanup/stats cleanup/dump cleanup/png


# This code removes straddler PRC/ERC voxels
function cleanup_prc()
{
  in=$1
  out=$2

  INFO=$(c3d $in -thresh 9 inf 1 0 -dilate 0 1x1x0vox  -o $TMPDIR/temp.nii.gz -info)
  NS=$(echo $INFO | sed -e "s/.*dim = .//g" -e "s/.;.*bb.*//g" | awk -F ',' '{print $3}')
  slicecmd=$(for((i=0;i<$NS;i++)); do echo "-push X -slice z $i -voxel-sum "; done)
  c3d $TMPDIR/temp.nii.gz -popas X $slicecmd | grep Voxel | awk '{print $3}' > $TMPDIR/counts.txt
  NNZ=$(cat $TMPDIR/counts.txt | grep -v '^0$' | wc -l)
  MEDIAN=$(cat $TMPDIR/counts.txt | grep -v '^0$' | sort -n | tail -n $((NNZ/2)) | head -n 1)
  CUTOFF=$((MEDIAN / 4))
  RULE=$(cat $TMPDIR/counts.txt | awk "{print NR-1,int(\$1 < $CUTOFF)}")
  c3d $TMPDIR/temp.nii.gz -cmv -replace $RULE -popas M $in -as X \
    -thresh 9 inf 1 0 -push M -times -scale -1 -shift 1 \
    -push X -times -o $out
  NLEFT=$(cat $TMPDIR/counts.txt | awk "\$1 > $CUTOFF {k++} END {print k}")
  echo $NLEFT
}

# Generate the statistics for the subject (similar to what's in the ASHS qc dir)
function genstats()
{
  in=$1
  c3d $in -dup -lstat > $TMPDIR/stat.txt
  THICK=$(c3d $in -info-full | grep Spacing | sed -e "s/[a-zA-Z:,]//g" -e "s/\]//" -e "s/\[//" | awk '{print $3}')
  VOLS=$(
  for i in 1 2 4 3 7 8 9 11 12 13; do
    cat $TMPDIR/stat.txt | awk "BEGIN {k=0;n=0} NR>1 && \$1 == $i { k=\$7; n=\$10; } END {print k,n}"
  done)
  echo $THICK,$VOLS | sed -e "s/ /,/g"
}

# Make PNG montage
function make_png()
{
  id=$1
  side=$2
  MRI=ashs/$id/tse_native_chunk_${side}.nii.gz
  SEG=cleanup/${id}_seg_${side}.nii.gz
  LABELS=/home/pauly/wolk/headtailatlas/snaplabels.txt

  # Generate coronal slices
  NSLICE=5
  for ((i=1; i<=$NSLICE; i++)); do

    PCT=$(echo $i $NSLICE | awk '{ print $1 * 100.0 / ($2 + 1) }')

    c3d $MRI -stretch 0 98% 0 255 -clip 0 255 -popas GG \
      $SEG -trim 40x40x0mm -as SS \
      -push GG -reslice-identity -push SS \
      -foreach -slice z ${PCT}% -flip xy -endfor \
      -popas S -popas G \
      -push G -type uchar -o $TMPDIR/cor_${id}_${side}_gray_${i}.png \
      -push S -oli $LABELS 0.5 -omc $TMPDIR/cor_${id}_${side}_seg_${i}.png

  done

  # Generate sagittal slices
  NSLICE=2
  for ((i=1; i<=$NSLICE; i++)); do

    PCT=$(echo $i $NSLICE | awk '{ print $1 * 100.0 / ($2 + 1) }')

    c3d $MRI -stretch 0 98% 0 255 -clip 0 255 -popas GG \
      $SEG -trim 0x40x40mm -resample 100x100x500% -as SS \
      -push GG -int 0 -reslice-identity -push SS \
      -foreach -slice x ${PCT}% -flip xy -endfor \
      -popas S -popas G \
      -push G -type uchar -o $TMPDIR/sag_${id}_${side}_gray_${i}.png \
      -push S -oli $LABELS 0.5 -omc $TMPDIR/sag_${id}_${side}_seg_${i}.png

  done

  montage \
    -tile 7x -geometry +5+5 \
    $TMPDIR/*_${id}_${side}_gray_*.png  $TMPDIR/*_${id}_${side}_seg_*.png \
    $TMPDIR/${id}_${side}_qa.png

  montage -label '%f' $TMPDIR/${id}_${side}_qa.png -geometry +1+1 \
    cleanup/png/${id}_${side}_qa.png
}

# Clean up an individual subject
function cleanup_subject()
{
  id=$1
  rm -rf cleanup/stats/stats_${id}.txt
  for side in left right; do
    cleanup_prc ashs/$id/bootstrap/fusion/lfseg_corr_usegray_${side}.nii.gz \
      cleanup/${id}_seg_${side}.nii.gz
    make_png $id $side
    echo $id,$side,$(genstats cleanup/${id}_seg_${side}.nii.gz) >> cleanup/stats/stats_${id}.txt
  done
}

# Clean up the names for the analysis
function cleanup_names()
{
  HEADER="ID,SEQ,SCANDATE,xval,ICV,side,Slice_Thickness"
  for sf in CA1 CA2 CA3 DG MISC SUB ERC BA35 BA36 CS; do
    HEADER="$HEADER,${sf}_vol,${sf}_ns"
  done
  echo $HEADER > stats_lr_cleanup.csv

  for fn in $(ls cleanup/stats/); do

    id=$(cat cleanup/stats/$fn | head -n 1 | awk -F ',' '{print $1}')

    # Get the ICV value
    ICV=$(cat ashs/${id}/final/${id}_icv.txt | awk '{print $2}')

    if [[ $(echo $id | grep xval) ]]; then

      dwid=$(echo $id | awk -F '_' '{print $3}')
      scandate=$(ls ~srdas/wd/ADC/$dwid/T00/rawNii \
        | grep 'DW.*nii.gz' | head -n 1 | awk -F '_' '{print $2}')

      newline="$dwid,T00,$scandate,xval,$ICV"

    else

      newline="$(echo $id | sed -e "s/_/,/g"),norm,$ICV"

    fi

    cat cleanup/stats/$fn | sed -e "s/$id/$newline/g" | tee -a stats_lr_cleanup.csv

  done
}

# MAIN
function main()
{
  for fn in $(cat cleanup_subj.txt); do  
    qsub -V -cwd -o cleanup/dump -j y -N "sfcleanup_${fn}" $0 cleanup_subject $fn
  done 

  qsub -V -cwd -o cleanup/dump -j y -hold_jid "sfcleanup_*" -sync y -b y sleep 1
  ### cleanup_names
}

if [[ ! $1 ]]; then
  main
elif [[ $1 = "cleanup_subject" ]]; then 
  cleanup_subject $2
fi
