#!/bin/bash 
#$ -S /bin/bash
set -x -e

# ----------
# ASHS Thickness Script
# ----------
#
# This is a thickness pipeline script for ASHS. It was used in preparation of the data
# for the 2014 Yushkevich et al. Human Brain Mapping paper 
# https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.22627
#
# This script is in rough shape. It works in our environment and it has been documented
# to help you run it in yours. But some modifications to the script may be necessary
# espeically if you use the script with a different ASHS atlas. This script is based on
# the ASHS PMC atlas used in the paper above.

# All comments that start with (!) are places that you may need to customize the script

# (!) Required input files for thickness analysis
# 
# analysis_input/subj.txt: a list of subject IDS separated by spaces, all on a single line
# 
# analysis_input/design_XXX.txt: design matrix for experiment XXX. You may have more than one
# experiment. Each row has the format
#   ID COVARIATE COVARIATE ... COVARIATE
# where ID is the subject id from subj.txt and covariates are things like Age, ICV, and group
# membership. For example, if you are doing a comparison between patients and controls and 
# covary for age and ICV your design might look like
#   Subj01 1 0 77 1002320
#   Subj02 0 1 75 991234
# here Subj01 is a patient and Subj02 is a control. For more info see documentation of meshglm
#
# analysis_input/contrast_XXX_YYY.txt: contrast vector for experiment XXX, contrast YYY. Each 
# experiment may have more than one contrast, although typically you would have one. For example
# you could for the same design matrix be interested in the contrast between patients and controls
# and regression with age. The contrast is a one-line file containing the contrast vector. For the
# age contrast, you would enter (for the design matrix above)
#   
#   0 0 1 0 
#
# and for the patient vs control contrast, you would enter
# 
#   1 -1 0 0 
#
# For more info about design matrices and contrasts, see docs for meshglm or SPM documentation.
#
# Note that XXX and YYY above can be arbitrary strings (without an underline)
#
# You also need to have processed a set of subjects using ASHS and to have Grid Engine (qsub)
# on your Linux machine. 

# (!) Location of ASHS. This must be the OLD (aka slow) ASHS that was in the 'master'
# branch of ASHS git as of 5/2018, not the new (aka fast) ASHS branch.
ASHS_ROOT=~/ashs

# (!) These are the ids of the labels for which thickness computation will be done.
# The labels combine some of the smaller labels, e.g., CA combines CA1, CA2 and CA3
LABEL_IDS=(BKG CA DG SUB ERC BA35 BA36 CS)

# (!) The list of labels in ASHS output corresponding to each of the labels above
LABEL_MRG=("0" "1 2 4" "3" "8" "9" "11" "12" "13")

# (!) The subset of labels that are relevant
LABEL_FG=(CA DG SUB ERC BA35 BA36)

# (!) This is the directory where you ran ASHS. The expectation is that for every subject
# there is a directory with that subject's ID that contains ASHS output. For each subject
# id, there should be a directory ${ASHSRUNDIR}/${id} containing the usual ASHS subdirs
# final, bootstrap, etc.
ASHSRUNDIR=~/shavg/parkinsons/ashs

# (!) This is the directory of ASHS segmentaion 'fix-up'. In our paper we cleaned up ASHS 
# segmentations by clearing slices that had very few segmented voxels. The fixup directory
# is expected to contain files with signature ${id}_seg_${side}.nii.gz where ${id} is the 
# subject id and ${side} is left or right
FIXUPDIR=~/exp04_headtail/fullset_truexval/cleanup

KINDS="tse mprage ${LABEL_IDS[*]}"
PATH=$ASHS_ROOT/ext/Linux/bin/ants_1042:$ASHS_ROOT/ext/Linux/bin:$PATH

mkdir -p dump

# This function copies subjects' data into the working directory
function copy_subject()
{
  fn=$1
  side=$2

  ### SEG=$ASHSRUNDIR/$fn/bootstrap/fusion/lfseg_heur_${side}.nii.gz
  SEG=$FIXUPDIR/${fn}_seg_${side}.nii.gz

  # Link the subfield images
  if [[ -f $SEG ]]; then
    ln -sf $ASHSRUNDIR/$fn/tse_to_chunktemp_${side}.nii.gz data/${fn}_${side}_tse.nii.gz
    ln -sf $ASHSRUNDIR/$fn/mprage_to_chunktemp_${side}.nii.gz data/${fn}_${side}_mprage.nii.gz
  fi

  # Generate a binary image for the label
  for ((i=0;i<${#LABEL_IDS[*]};i++)); do

    c3d $SEG -replace $(for k in ${LABEL_MRG[i]}; do echo $k 999; done) -thresh 999 999 1 0 \
      -o $TMPDIR/binary_${LABEL_IDS[i]}.nii.gz

    WarpImageMultiTransform 3  $TMPDIR/binary_${LABEL_IDS[i]}.nii.gz data/${fn}_${side}_${LABEL_IDS[i]}.nii.gz \
      -R data/${fn}_${side}_tse.nii.gz \
      $ASHSRUNDIR/$fn/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz \
      $ASHSRUNDIR/$fn/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
      $ASHSRUNDIR/$fn/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

  done

  # Vote in the subject space
  c3d $(for sub in ${LABEL_IDS[*]}; do echo data/${fn}_${side}_${sub}.nii.gz; done) \
    -vote -type ushort -o data/${fn}_${side}_seg.nii.gz
}

# This function runs qsub to copy all subjects' data
function copy_data() 
{
  # Get the data
  mkdir -p data
  for id in $(cat analysis_input/subj.txt); do
    fn=$(ls ashs | grep $id)
    for side in left right; do
    
      qsub -V -cwd -o dump -j y -N "thickinit_${fn}_${side}" $0 copy_subject $fn $side

    done
  done

  qsub -V -cwd -o dump -j y -hold_jid "thickinit_*" -sync y -b y sleep 1
}


function average_subfield()
{
  side=$1
  kind=$2

  AverageImages 3 work/template_fullchunk_${side}_${kind}.nii.gz 0 data/*_${side}_${kind}.nii.gz
}


function initial_average()
{
  mkdir -p work


  # Compute initial average for each subfield mask
  for side in left right; do

    # Average all the input images
    for kind in $KINDS; do
      average_subfield $side $kind
    done

    # Compute the initial mask for the segmentation
    c3d $(for sub in ${LABEL_IDS[*]}; do echo work/template_fullchunk_${side}_${sub}.nii.gz; done | grep -v BKG) \
      -mean -thresh 1e-5 inf 1 0 -trim 5vox \
      -o work/template_mask_${side}.nii.gz

    # Trim every template component using the mask
    for kind in $KINDS; do

      # Old code (use the average as the template)
      c3d work/template_mask_${side}.nii.gz \
        work/template_fullchunk_${side}_${kind}.nii.gz \
        -reslice-identity -o work/template_${side}_${kind}.nii.gz

      # (!) if you want to use a specific subject as the template rather than the shape average of all
      # subjects, uncomment the lines below. INIT should point to a subject ID
      ### if [[ side=='left' ]]; then INIT=DW228; else INIT=DW255; fi
      ### id=$(ls $ASHSRUNDIR | grep $INIT)
      ### c3d work/template_mask_${side}.nii.gz \
      ###   data/${id}_${side}_${kind}.nii.gz \
      ###   -reslice-identity -o work/template_${side}_${kind}.nii.gz

    done

    

  done


}

function shape_update_to_template()
{
  side=$1

  # Borrowed from ANTS buildtemplateparallel

  # Average the warp fields
  local TEMPWARP=work/template_${side}_warp.nii.gz
  AverageImages 3 $TEMPWARP 0 work/*_${side}_totempWarp.nii.gz

  # Scale by -0.25 (gradient step)
  MultiplyImages 3 $TEMPWARP -0.25 $TEMPWARP

  # Create a new template affine
  local TEMPAFF=work/template_${side}_Affine.txt
  if [[ -f $TEMPAFF ]]; then rm -rf $TEMPAFF; fi

  cat work/*_${side}_totempAffine.txt | grep '^Parameters:' | awk '\
    BEGIN { for(i=0;i<12;i++) x[i]=0} \
    { for(i=0;i<12;i++) x[i]=x[i]+$(i+2) } \
    END { \
      printf "Transform: MatrixOffsetTransformBase_double_3_3\nParameters: "; \
      for(i=0;i<12;i++) printf "%f ",x[i]/NR; \
      printf "\nFixedParameters: 0 0 0\n";}' > $TEMPAFF

  # Compose the warps
  WarpImageMultiTransform 3 $TEMPWARP $TEMPWARP -i $TEMPAFF -R work/template_${side}_tse.nii.gz

  TEMPWARPFULL=work/template_${side}_fullwarp.nii.gz
  ComposeMultiTransform 3 \
    $TEMPWARPFULL -R work/template_${side}_tse.nii.gz \
    -i $TEMPAFF $TEMPWARP $TEMPWARP $TEMPWARP $TEMPWARP

  # Apply this warp to all the template derivatives
  for kind in $KINDS; do

    WarpImageMultiTransform 3 \
      work/template_${side}_${kind}.nii.gz \
      work/template_${side}_${kind}.nii.gz \
      $TEMPWARPFULL

  done
}
      
function ants_iter()
{
  id=$1
  side=$2
  doants=$3

  # Before we vote, use ml_affine for nice affine alignment
  ml_affine \
    work/template_${side}_seg.nii.gz \
    data/${id}_${side}_seg.nii.gz \
    work/${id}_${side}_mlaffine.txt

  # Do we do ANTS or not?
  if [[ $doants -eq 0 ]]; then

    for sub in $KINDS; do

      c3d \
        work/template_${side}_tse.nii.gz \
        data/${id}_${side}_${sub}.nii.gz \
        -reslice-matrix work/${id}_${side}_mlaffine.txt \
        -o work/${id}_${side}_totemp_reslice_${sub}.nii.gz 

    done

  else

    # Convert that to ITK format
    c3d_affine_tool work/${id}_${side}_mlaffine.txt -oitk work/${id}_${side}_mlaffine_itk.txt

    local FIX=work/template_${side}_seg.nii.gz
    local MOV=data/${id}_${side}_seg.nii.gz

    # Long's parameters
    WGT=$(echo ${#LABEL_IDS[*]} | awk '{print 1.0 / $1}')
    CMD=""
    for sub in ${LABEL_IDS[*]}; do
      CMD="$CMD -m MSQ[work/template_${side}_${sub}.nii.gz,data/${id}_${side}_${sub}.nii.gz,1]"
    done

    ANTS 3 \
      $CMD \
      -t SyN[0.25] -r Gauss[0.5,0] -i 80x80x20 \
      -x work/template_mask_${side}.nii.gz \
      -a work/${id}_${side}_mlaffine_itk.txt \
      --continue-affine 0 \
      --use-all-metrics-for-convergence \
      -o work/${id}_${side}_totemp.nii.gz | tee work/${id}_${side}_antsoutput.txt
    
    for sub in $KINDS; do

      WarpImageMultiTransform 3 \
        data/${id}_${side}_${sub}.nii.gz \
        work/${id}_${side}_totemp_reslice_${sub}.nii.gz \
        -R work/template_${side}_tse.nii.gz \
        work/${id}_${side}_totempWarp.nii.gz \
        work/${id}_${side}_totempAffine.txt

    done

    c3d $(for sub in ${LABEL_IDS[*]}; do echo work/${id}_${side}_totemp_reslice_${sub}.nii.gz; done) \
      -vote -type ushort -o work/${id}_${side}_totemp_reslice_seg.nii.gz

  fi
}

function main_loop
{
  # Main iteration loop
  NITER=8
  SITER=0
  for side in left right; do

    IDS=$(ls data | grep ${side}_tse | sed -e "s/_${side}_tse.nii.gz//")

    for ((iter=$SITER;iter<$NITER;iter++)); do

      # If not the first iteration, average the resliced posterior maps

      # Create the segmentation for the template
      c3d $(for sub in ${LABEL_IDS[*]}; do echo work/template_${side}_${sub}.nii.gz; done) \
        -vote -type ushort -o work/template_${side}_seg.nii.gz

      # Store the segmentation for posterity
      cp -a work/template_${side}_seg.nii.gz work/template_${side}_iter${iter}_seg.nii.gz

      # Back up template
      ITDIR=work/$(printf iter_%s_%02d $side $iter)
      mkdir -p $ITDIR
      cp -a work/template_${side}_*.nii.gz $ITDIR/

      # Do ants?
      if [[ iter -lt 3 ]]; then doants=0; else doants=1; fi

      # Run ANTS for each image
      for id in $IDS; do

        # Submit ANTS job
        qsub -V -cwd -o dump -j y -N "tseants_${id}_${side}_iter${iter}" $0 ants_iter $id $side $doants

      done

      # Wait for completion
      qsub -V -cwd -o dump -j y -hold_jid "tseants*" -sync y -b y sleep 1

      # If this is the last iteration, we don't want to recompute the template
      if [[ $iter -lt $((NITER-1)) ]]; then

        # Compute average images
        for kind in $KINDS; do

          if [[ $kind = "tse" ]]; then NORM=1; else NORM=0; fi
          AverageImages 3 work/template_${side}_${kind}.nii.gz \
            $NORM work/*_${side}_totemp_reslice_${kind}.nii.gz

        done

        # Perform shape averaging
        if [[ $doants -eq 1 ]]; then
          shape_update_to_template $side
        fi
      fi

    done
  done
}



# Extract meshes for the subfields and apply warps to these meshes. This allows us to perform
# statistical analysis on the mesh boundaries
function warp_meshes()
{
  mkdir -p meshwarp
  mkdir -p jacobian

  # Iterate over side and subject
  for side in left right; do

    # Generate meshes for the individual subfields
    for sub in ${LABEL_FG[*]}; do

      vtklevelset work/template_${side}_${sub}.nii.gz meshwarp/template_${side}_${sub}.vtk 0.5

    done

    # (!) Generate one compound label for all the non-DG subfields. This code is specific to the 
    # set of labels used in the ASHS PMC atlas, and generates a single 'MRG' label that combines
    # CA, SUB, ERC, BA35, BA36 labels. If you are using a different set of subfields, change this
    # line accordingly.
    c3d \
      work/template_${side}_CA.nii.gz \
      work/template_${side}_SUB.nii.gz \
      work/template_${side}_ERC.nii.gz \
      work/template_${side}_BA35.nii.gz \
      work/template_${side}_BA36.nii.gz \
      work/template_${side}_BKG.nii.gz -scale -1 \
      work/template_${side}_CS -scale -1 \
      work/template_${side}_DG.nii.gz -scale -1 \
      -add -add -add -add -add -add -add -o /tmp/mergesf.nii.gz

    vtklevelset /tmp/mergesf.nii.gz meshwarp/template_${side}_MRG.vtk 0.0

     
    IDS=$(ls data | grep ${side}_tse | sed -e "s/_${side}_tse.nii.gz//")
    
    for id in $IDS; do

      # Submit job for this subject
      qsub -V -cwd -o dump -j y -N "warpmesh_${id}" $0 warp_meshes_subj $id $side

    done

  done
  
  # Wait for completion
  qsub -V -cwd -o dump -j y -hold_jid "warpmesh*" -sync y -b y sleep 1

}

function warp_meshes_subj()
{
  id=$1
  side=$2
  ALLSF="${LABEL_FG[*]} MRG"

  # Generate and save the full warps

  # Compose the transformation between the template and the subject
  ComposeMultiTransform 3 \
    $TMPDIR/compose.nii \
    -R work/template_${side}_tse.nii.gz \
    work/${id}_${side}_totempWarp.nii.gz \
    work/${id}_${side}_totempAffine.txt \
    ashs/${id}/ants_t1_to_temp/ants_t1_to_tempWarp.nii \
    ashs/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
    ashs/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt


  # Split the transformations into xyz
  c3d -mcs $TMPDIR/compose.nii -oo $TMPDIR/comp%02d.nii

  # Apply the warp to each of the meshes
  for sub in $ALLSF; do

    # Warp the subfield into subject space
    warpmesh -w ants -m ras \
      meshwarp/template_${side}_${sub}.vtk \
      meshwarp/${id}_${side}_${sub}_tempfit.vtk \
      $TMPDIR/comp00.nii $TMPDIR/comp01.nii $TMPDIR/comp02.nii

    # Extract the thickness of the subfield
    cmrep_vskel -Q ~/bin/qvoronoi \
      -T meshwarp/${id}_${side}_${sub}_thickmap.vtk -p 1.2 -e 6 \
      meshwarp/${id}_${side}_${sub}_tempfit.vtk $TMPDIR/skel.vtk 

  done

  # Compute the Jacobian map
  ### ANTSJacobian 3 $TMPDIR/compose.nii jacobian/${id}_${side}_totemp 1
}

# Warp template segmentations to subject space
function warp_tosubj()
{
  mkdir -p tosubj

  # Iterate over side and subject
  for side in left right; do

    IDS=$(ls data | grep ${side}_tse | sed -e "s/_${side}_tse.nii.gz//")
    
    for id in $IDS; do

      # Submit job for this subject
      qsub -V -cwd -o dump -j y -N "warpsubj_${id}" $0 warp_tosubj_subj $id $side

    done

  done
  
  # Wait for completion
  qsub -V -cwd -o dump -j y -hold_jid "warpsubj*" -sync y -b y sleep 1

}

function warp_tosubj_subj()
{
  id=$1
  side=$2

  # We also want to compute smoothed subfield masks and volumes
  c3d ashs/$id/tse_native_chunk_${side}.nii.gz -resample 100x100x500% \
    -type short -o $TMPDIR/refspace.nii.gz

  ComposeMultiTransform 3 \
    $TMPDIR/compinv.nii \
    -R $TMPDIR/refspace.nii.gz \
    -i ashs/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
    -i ashs/${id}/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
    ashs/${id}/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii \
    -i work/${id}_${side}_totempAffine.txt \
    work/${id}_${side}_totempInverseWarp.nii.gz

  for ((i=0;i<${#LABEL_IDS[*]};i++)); do

    sub=${LABEL_IDS[i]}

    WarpImageMultiTransform 3 \
      work/template_${side}_${sub}.nii.gz \
      $TMPDIR/reslice_${sub}.nii.gz \
      -R $TMPDIR/refspace.nii.gz \
      $TMPDIR/compinv.nii 

  done

  mkdir -p tosubj
  c3d $(for sub in ${LABEL_IDS[*]}; do echo $TMPDIR/reslice_${sub}.nii.gz; done) \
    -vote -type ushort -o tosubj/${id}_${side}_seg.nii.gz
}

function disp_stats()
{
  ALLSF="${LABEL_FG[*]} MRG"
  mkdir -p meshwarp/analysis
  for side in left right; do
    for sub in $ALLSF; do

      # Submit job for this subject
      qsub -V -cwd -o dump -j y -N "meshstat_${side}_${sub}" $0 disp_stats_qsub $side $sub

    done
  done
  
  # Wait for completion
  qsub -V -cwd -o dump -j y -hold_jid "meshstat_*" -sync y -b y sleep 1
}

function disp_stats_qsub()
{
  side=$1
  sub=$2

  # Displacement analysis
  MESHES=$(for id in $(echo $(cat analysis_input/subj.txt)); do \
    echo $(find meshwarp | grep $id | grep $side | grep $sub | grep tempfit); done)

  meshdisp $MESHES meshwarp/analysis/disp_${side}_${sub}.vtk

  # Thickness analysis
  MESHES=$(for id in $(echo $(cat analysis_input/subj.txt)); do \
    echo $(find meshwarp | grep $id | grep $side | grep $sub | grep thickmap); done)

  mesh_merge_arrays -r meshwarp/template_${side}_${sub}.vtk \
    meshwarp/analysis/thick_${side}_${sub}.vtk Thickness $MESHES

}

function average_segmentation_error_map()
{
  mkdir -p segerrormap

  for side in left right; do
    # Compute an average disagreement map
    for fn in $(ls ashs | grep xval); do

      # Compute a difference map
      c3d \
        ashs/$fn/bootstrap/fusion/lfseg_corr_usegray_${side}.nii.gz -dup \
        ashs/$fn/refseg/refseg_${side}.nii.gz -int 0 -reslice-identity \
        -scale -1 -add -thresh 0 0 0 1 -o /tmp/diff.nii.gz

      # Warp the difference map into the template space
      WarpImageMultiTransform 3 \
        /tmp/diff.nii.gz /tmp/diffmap_${fn}_${side}.nii.gz -R work/template_${side}_tse.nii.gz \
        ashs/$fn/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz \
        ashs/$fn/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
        ashs/$fn/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

    done

    AverageImages 3 segerrormap/segerrormap_${side}.nii.gz 0 /tmp/diffmap_*_${side}.nii.gz
  done
}

# (!) For thickness analysis, we can run analysis on multiple groups of meshes at once. For example
# in the first experiment, we merge all the non-DG structures into one strip. In the second, we 
# separately analyze the thickness of different subfields
ANGRP[0]="MRG DG"
ANGRP[1]="CA DG SUB ERC BA35 BA36"

GRPNM[0]="merged"
GRPNM[1]="all"

function thick_stats()
{
  for design in $(ls analysis_input | grep "design_.*txt"); do

    exp=$(echo $design | sed -e "s/^design_//" -e "s/\.txt//")

    for ((igrp=0;igrp<${#ANGRP[*]};igrp++)); do

      # Submit job for this subject
      qsub -V -cwd -o dump -j y -N "thickstat_${exp}_${GRPNM[igrp]}" $0 thick_stats_qsub $exp $igrp

    done

  done
  
  # Wait for completion
  qsub -V -cwd -o dump -j y -hold_jid "thickstat_*" -sync y -b y sleep 1
}

function thick_stats_qsub()
{
  exp=$1
  igrp=$2

  MYGRP=${ANGRP[igrp]}
  GNAME=${GRPNM[igrp]}

  # Create the work directory for this analysis
  WORK=meshwarp/analysis/design_${exp}_group_${GNAME}
  mkdir -p $WORK
  rm -rf $WORK/*

  # Get the list of subjects
  SUBJ=$(cat analysis_input/design_${exp}.txt | awk '{print $1}')

  # Merge the meshes for this analysis
  for side in left right; do
    for sub in $MYGRP; do

      MESHES=$(for id in $SUBJ; do \
        echo $(find meshwarp | grep "${id}.*_${side}_${sub}_thickmap.vtk"); done)

      mesh_merge_arrays -r meshwarp/template_${side}_${sub}.vtk \
        $WORK/thick_${side}_${sub}.vtk Thickness $MESHES

    done
  done

  # Generate the design matrix for meshglm
  cat analysis_input/design_${exp}.txt | awk '{$1=""; print}' > $WORK/design.txt

  # Go through the list of contrasts
  for con in $(ls analysis_input | grep "contrast_${exp}_.*\.txt"); do

    # Get the suffix for this contrast
    suffix=$(echo $con | sed -e "s/^contrast_${exp}_//" -e "s/\.txt//")

    # Create the directory for this contrast
    CWORK=$WORK/contrast_${suffix}
    mkdir -p $CWORK
    
    FULLNM="design_${exp}_group_${GNAME}_con_${suffix}"

    # Copy the contrast
    cp analysis_input/$con $CWORK/contrast.txt

    # Build the list of meshes to include
    MESHPARAM=""
    for side in left right; do
      for sub in $MYGRP; do
        MESHPARAM="$MESHPARAM -m $WORK/thick_${side}_${sub}.vtk $CWORK/thickstat_${FULLNM}_${side}_${sub}.vtk"

      done
    done

    meshglm $MESHPARAM \
      -g $WORK/design.txt $CWORK/contrast.txt \
      -a Thickness -d 10

  done
}

function reset_dir()
{
  rm -rf data work meshwarp jacobian segerrormap tosubj 
  rm -rf dump/*
}


if [[ $1 == "copy_subject" ]]; then

  copy_subject $2 $3

elif [[ $1 == "average_subfield" ]]; then

  average_subfield $2 $3

elif [[ $1 == "ants_iter" ]]; then

  ants_iter $2 $3 $4

elif [[ $1 == "warp_meshes_subj" ]]; then

  warp_meshes_subj $2 $3

elif [[ $1 == "warp_tosubj_subj" ]]; then

  warp_tosubj_subj $2 $3

elif [[ $1 == "disp_stats_qsub" ]]; then

  disp_stats_qsub $2 $3

elif [[ $1 == "thick_stats_qsub" ]]; then

  thick_stats_qsub $2 $3

else

  # (!) Uncomment this line to clean up the output directories
  ### reset_dir
  copy_data
  initial_average
  main_loop
  warp_meshes
  warp_tosubj
  disp_stats
  thick_stats

  # (!) this is optional, computes a map of segmentation error
  ### #average_segmentation_error_map
  exit

fi

