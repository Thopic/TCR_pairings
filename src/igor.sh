N=100000000

# BATCHNAME=beta_human_werr_1e8
# CMD="./../../../IGoR/bin/igor -set_wd ../Datas/IGoR_gen/ -batch $BATCHNAME"
# $CMD -species human -chain beta -generate $N --CDR3
# # $CMD -read_seqs "../Datas/IGoR_gen/${BATCHNAME}_generated/generated_seqs_noerr.csv"
# # $CMD -species human -chain beta -align --all
# # $CMD -species human -chain beta -evaluate -output --Pgen

BATCHNAME=alpha_human_werr_1e8
CMD="./../../../IGoR/bin/igor -set_wd ../Datas/IGoR_gen/ -batch $BATCHNAME"
$CMD -species human -chain alpha -generate $N --CDR3
# $CMD -read_seqs "../Datas/IGoR_gen/${BATCHNAME}_generated/generated_seqs_noerr.csv"
# $CMD -species human -chain alpha -align --all
# $CMD -species human -chain alpha -evaluate -output --Pgen
