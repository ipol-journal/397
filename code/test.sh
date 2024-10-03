imname=$1

input_image=../images/$imname.png
output_dir=../images/results/$imname

if [ ! -d "$output_dir" ]
then
    mkdir $output_dir
fi

########################################################################
# y
###############

if [ "$imname" = "y" ]
then
    operator=opening #erosion
    sigma=5
    p=7
    ang_prec=9 # in [0, 90] degrees
    cone_ang=30 # in [0, 90] degrees
    n_it=10
    apply_to=anisotropy # anisotropy or input
    mode=bin # bin or nbin
fi

########################################################################

########################################################################
# k
###############

if [ "$imname" = "k" ]
then
    operator=opening
    sigma=5
    p=7
    ang_prec=9 # in [0, 90] degrees
    cone_ang=30 # in [0, 90] degrees
    n_it=15
    apply_to=anisotropy # anisotropy or input
    mode=bin # bin or nbin
fi

########################################################################

########################################################################
# retine
###############

if [ "$imname" = "retine" ]
then
    operator=opening
    sigma=2
    p=3
    ang_prec=15 # in [0, 90] degrees
    cone_ang=30 # in [0, 90] degrees
    n_it=5
    apply_to=anisotropy # anisotropy or input
    mode=bin # bin or nbin
fi

########################################################################

########################################################################
# eye crop
###############

if [ "$imname" = "eye_crop" ]
then
    operator=asfclos
    sigma=5
    p=7
    ang_prec=9 # in [0, 90] degrees
    cone_ang=30 # in [0, 90] degrees
    n_it=3
    apply_to=input # anisotropy or input
    mode=bin # bin or nbin
fi

########################################################################

########################################################################
# autumn
###############

if [ "$imname" = "autumn" ]
then
    operator=asfclos
    sigma=2
    p=3
    ang_prec=9 # in [0, 90] degrees
    cone_ang=30 # in [0, 90] degrees
    n_it=5
    apply_to=input # anisotropy or input
    mode=nbin # bin or nbin
fi

########################################################################

########################################################################
# klok
###############

if [ "$imname" = "klok" ]
then
    operator=closing
    sigma=4
    p=7
    ang_prec=9 # in [0, 90] degrees
    cone_ang=30 # in [0, 90] degrees
    n_it=2
    apply_to=input # anisotropy or input
    mode=bin # bin or nbin
fi

########################################################################

########################################################################
# klok
###############

if [ "$imname" = "klok" ]
then
    operator=closing
    sigma=4
    p=7
    ang_prec=9 # in [0, 90] degrees
    cone_ang=30 # in [0, 90] degrees
    n_it=2
    apply_to=input # anisotropy or input
    mode=bin # bin or nbin
fi

########################################################################

output_image=$output_dir/$imname\_$operator\_s$sigma\_p$p\_alpha$ang_prec\_beta$cone_ang\_it$n_it\_$mode.png

anisotropy_image=$output_dir/$imname\_anisotropy_s$sigma.png

./build/anisop $operator $input_image $output_image $anisotropy_image $sigma $p $ang_prec $cone_ang $n_it $apply_to $mode


