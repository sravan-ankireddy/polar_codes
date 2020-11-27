#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "functions_SCL.h"

int main(){

    // srand(time(NULL));

	clock_t start, end;
    clock_t enc_start, enc_end;
    clock_t dec_start_rec, dec_end_rec;
    clock_t dec_start, dec_end;

    float cpu_time_used;
    float enc_cpu_time_used = 0;
    float dec_cpu_time_used_rec = 0;
    float dec_cpu_time_used = 0;
    start = clock();

/* Channel reliability in increasing order*/
int Q[1024] = {0 ,1 ,2 ,4 ,8 ,16 ,32 ,3 ,5 ,64 ,9 ,6 ,17 ,10 ,18 ,128 ,12 ,33 ,65 ,20 ,256 ,34 ,24 ,36 ,7 ,129 ,66 ,512 ,11 ,40 ,68 ,130 ,
		19 ,13 ,48 ,14 ,72 ,257 ,21 ,132 ,35 ,258 ,26 ,513 ,80 ,37 ,25 ,22 ,136 ,260 ,264 ,38 ,514 ,96 ,67 ,41 ,144 ,28 ,69 ,42 ,
		516 ,49 ,74 ,272 ,160 ,520 ,288 ,528 ,192 ,544 ,70 ,44 ,131 ,81 ,50 ,73 ,15 ,320 ,133 ,52 ,23 ,134 ,384 ,76 ,137 ,82 ,56 ,27 ,
		97 ,39 ,259 ,84 ,138 ,145 ,261 ,29 ,43 ,98 ,515 ,88 ,140 ,30 ,146 ,71 ,262 ,265 ,161 ,576 ,45 ,100 ,640 ,51 ,148 ,46 ,75 ,266 ,273 ,517 ,104 ,162 ,
		53 ,193 ,152 ,77 ,164 ,768 ,268 ,274 ,518 ,54 ,83 ,57 ,521 ,112 ,135 ,78 ,289 ,194 ,85 ,276 ,522 ,58 ,168 ,139 ,99 ,86 ,60 ,280 ,89 ,290 ,529 ,524 ,
		196 ,141 ,101 ,147 ,176 ,142 ,530 ,321 ,31 ,200 ,90 ,545 ,292 ,322 ,532 ,263 ,149 ,102 ,105 ,304 ,296 ,163 ,92 ,47 ,267 ,385 ,546 ,324 ,208 ,386 ,150 ,153 ,
		165 ,106 ,55 ,328 ,536 ,577 ,548 ,113 ,154 ,79 ,269 ,108 ,578 ,224 ,166 ,519 ,552 ,195 ,270 ,641 ,523 ,275 ,580 ,291 ,59 ,169 ,560 ,114 ,277 ,156 ,87 ,197 ,
		116 ,170 ,61 ,531 ,525 ,642 ,281 ,278 ,526 ,177 ,293 ,388 ,91 ,584 ,769 ,198 ,172 ,120 ,201 ,336 ,62 ,282 ,143 ,103 ,178 ,294 ,93 ,644 ,202 ,592 ,323 ,392 ,
		297 ,770 ,107 ,180 ,151 ,209 ,284 ,648 ,94 ,204 ,298 ,400 ,608 ,352 ,325 ,533 ,155 ,210 ,305 ,547 ,300 ,109 ,184 ,534 ,537 ,115 ,167 ,225 ,326 ,306 ,772 ,157 ,
		656 ,329 ,110 ,117 ,212 ,171 ,776 ,330 ,226 ,549 ,538 ,387 ,308 ,216 ,416 ,271 ,279 ,158 ,337 ,550 ,672 ,118 ,332 ,579 ,540 ,389 ,173 ,121 ,553 ,199 ,784 ,179 ,
		228 ,338 ,312 ,704 ,390 ,174 ,554 ,581 ,393 ,283 ,122 ,448 ,353 ,561 ,203 ,63 ,340 ,394 ,527 ,582 ,556 ,181 ,295 ,285 ,232 ,124 ,205 ,182 ,643 ,562 ,286 ,585 ,
		299 ,354 ,211 ,401 ,185 ,396 ,344 ,586 ,645 ,593 ,535 ,240 ,206 ,95 ,327 ,564 ,800 ,402 ,356 ,307 ,301 ,417 ,213 ,568 ,832 ,588 ,186 ,646 ,404 ,227 ,896 ,594 ,
		418 ,302 ,649 ,771 ,360 ,539 ,111 ,331 ,214 ,309 ,188 ,449 ,217 ,408 ,609 ,596 ,551 ,650 ,229 ,159 ,420 ,310 ,541 ,773 ,610 ,657 ,333 ,119 ,600 ,339 ,218 ,368 ,
		652 ,230 ,391 ,313 ,450 ,542 ,334 ,233 ,555 ,774 ,175 ,123 ,658 ,612 ,341 ,777 ,220 ,314 ,424 ,395 ,673 ,583 ,355 ,287 ,183 ,234 ,125 ,557 ,660 ,616 ,342 ,316 ,
		241 ,778 ,563 ,345 ,452 ,397 ,403 ,207 ,674 ,558 ,785 ,432 ,357 ,187 ,236 ,664 ,624 ,587 ,780 ,705 ,126 ,242 ,565 ,398 ,346 ,456 ,358 ,405 ,303 ,569 ,244 ,595 ,
		189 ,566 ,676 ,361 ,706 ,589 ,215 ,786 ,647 ,348 ,419 ,406 ,464 ,680 ,801 ,362 ,590 ,409 ,570 ,788 ,597 ,572 ,219 ,311 ,708 ,598 ,601 ,651 ,421 ,792 ,802 ,611 ,
		602 ,410 ,231 ,688 ,653 ,248 ,369 ,190 ,364 ,654 ,659 ,335 ,480 ,315 ,221 ,370 ,613 ,422 ,425 ,451 ,614 ,543 ,235 ,412 ,343 ,372 ,775 ,317 ,222 ,426 ,453 ,237 ,
		559 ,833 ,804 ,712 ,834 ,661 ,808 ,779 ,617 ,604 ,433 ,720 ,816 ,836 ,347 ,897 ,243 ,662 ,454 ,318 ,675 ,618 ,898 ,781 ,376 ,428 ,665 ,736 ,567 ,840 ,625 ,238 ,
		359 ,457 ,399 ,787 ,591 ,678 ,434 ,677 ,349 ,245 ,458 ,666 ,620 ,363 ,127 ,191 ,782 ,407 ,436 ,626 ,571 ,465 ,681 ,246 ,707 ,350 ,599 ,668 ,790 ,460 ,249 ,682 ,
		573 ,411 ,803 ,789 ,709 ,365 ,440 ,628 ,689 ,374 ,423 ,466 ,793 ,250 ,371 ,481 ,574 ,413 ,603 ,366 ,468 ,655 ,900 ,805 ,615 ,684 ,710 ,429 ,794 ,252 ,373 ,605 ,
		848 ,690 ,713 ,632 ,482 ,806 ,427 ,904 ,414 ,223 ,663 ,692 ,835 ,619 ,472 ,455 ,796 ,809 ,714 ,721 ,837 ,716 ,864 ,810 ,606 ,912 ,722 ,696 ,377 ,435 ,817 ,319 ,
		621 ,812 ,484 ,430 ,838 ,667 ,488 ,239 ,378 ,459 ,622 ,627 ,437 ,380 ,818 ,461 ,496 ,669 ,679 ,724 ,841 ,629 ,351 ,467 ,438 ,737 ,251 ,462 ,442 ,441 ,469 ,247 ,
		683 ,842 ,738 ,899 ,670 ,783 ,849 ,820 ,728 ,928 ,791 ,367 ,901 ,630 ,685 ,844 ,633 ,711 ,253 ,691 ,824 ,902 ,686 ,740 ,850 ,375 ,444 ,470 ,483 ,415 ,485 ,905 ,
		795 ,473 ,634 ,744 ,852 ,960 ,865 ,693 ,797 ,906 ,715 ,807 ,474 ,636 ,694 ,254 ,717 ,575 ,913 ,798 ,811 ,379 ,697 ,431 ,607 ,489 ,866 ,723 ,486 ,908 ,718 ,813 ,
		476 ,856 ,839 ,725 ,698 ,914 ,752 ,868 ,819 ,814 ,439 ,929 ,490 ,623 ,671 ,739 ,916 ,463 ,843 ,381 ,497 ,930 ,821 ,726 ,961 ,872 ,492 ,631 ,729 ,700 ,443 ,741 ,
		845 ,920 ,382 ,822 ,851 ,730 ,498 ,880 ,742 ,445 ,471 ,635 ,932 ,687 ,903 ,825 ,500 ,846 ,745 ,826 ,732 ,446 ,962 ,936 ,475 ,853 ,867 ,637 ,907 ,487 ,695 ,746 ,
		828 ,753 ,854 ,857 ,504 ,799 ,255 ,964 ,909 ,719 ,477 ,915 ,638 ,748 ,944 ,869 ,491 ,699 ,754 ,858 ,478 ,968 ,383 ,910 ,815 ,976 ,870 ,917 ,727 ,493 ,873 ,701 ,
		931 ,756 ,860 ,499 ,731 ,823 ,922 ,874 ,918 ,502 ,933 ,743 ,760 ,881 ,494 ,702 ,921 ,501 ,876 ,847 ,992 ,447 ,733 ,827 ,934 ,882 ,937 ,963 ,747 ,505 ,855 ,924 ,
		734 ,829 ,965 ,938 ,884 ,506 ,749 ,945 ,966 ,755 ,859 ,940 ,830 ,911 ,871 ,639 ,888 ,479 ,946 ,750 ,969 ,508 ,861 ,757 ,970 ,919 ,875 ,862 ,758 ,948 ,977 ,923 ,
		972 ,761 ,877 ,952 ,495 ,703 ,935 ,978 ,883 ,762 ,503 ,925 ,878 ,735 ,993 ,885 ,939 ,994 ,980 ,926 ,764 ,941 ,967 ,886 ,831 ,947 ,507 ,889 ,984 ,751 ,942 ,996 ,
		971 ,890 ,509 ,949 ,973 ,1000 ,892 ,950 ,863 ,759 ,1008 ,510 ,979 ,953 ,763 ,974 ,954 ,879 ,981 ,982 ,927 ,995 ,765 ,956 ,887 ,985 ,997 ,986 ,943 ,891 ,998 ,766 ,
		511 ,988 ,1001 ,951 ,1002 ,893 ,975 ,894 ,1009 ,955 ,1004 ,1010 ,957 ,983 ,958 ,987 ,1012 ,999 ,1016 ,767 ,989 ,1003 ,990 ,1005 ,959 ,1011 ,1013 ,895 ,1006 ,1014 ,1017 ,1018 ,
		991 ,1020 ,1007 ,1015 ,1019 ,1021 ,1022 ,1023};

/* Code Parameters */

	/* Length of code */
	int N = 1024;

    /* Length of CRC */
    int crc_l = 8;

    /* CRC polynomial */
    int polynomial[9] = {1, 1, 1, 0, 1, 0, 1, 0, 1};

    /* Depth of tree */
    int n = 0, N_temp = N;
    
    /* log function */
    while (N_temp >>= 1) n++;

	/* Rate of code */
	float rate = 0.5;

	/* List size -- 1,2,4,8,16 */
	int l = 16;

	/* Number of information bits */
	int K = (int)N*rate;

    /* Boolean array with information nodes pos = 1 */
    int info_nodes[N];

	/* Position of frozen bits */
	int frozen_pos[N-K];
	int i_Q = 0;
	for (i_Q = 0; i_Q < N-K; i_Q++)
	{
		frozen_pos[i_Q] = Q[i_Q];
        info_nodes[Q[i_Q]] = 0;
	}

	/* Position of Information bits */
	int data_pos[K];

	for (i_Q = 0; i_Q < K; i_Q++)
	{
		data_pos[i_Q] = Q[i_Q + N-K];
        info_nodes[Q[i_Q + N-K]] = 1;
	}

/* Simulation Parameters */
	
	/* Number of Simulations */
	int num_sim = 100*1000;

	/* No. of levels of Noise */
	int num_EbN0dB = 10;

    /* Max received value */
    int rmax = 3;

    /* Max integer received values */
    int maxqr = 31;
	
	/* Eb/N in dB */
	float EbN0dB[num_EbN0dB];
	int i_e = 0;
    EbN0dB[0] = 1;
    for (i_e = 1; i_e < num_EbN0dB; i_e++)
    {
        EbN0dB[i_e] = EbN0dB[i_e - 1] + 0.3;
    }

    /* Standard Deviation of noise */
	float sig;		
	float sigma[num_EbN0dB];
    int i_std = 0;
    for (i_std = 0; i_std<num_EbN0dB; i_std++)
    {
        sigma[i_std] = sqrt(1/(2 * rate) * pow(10,(-EbN0dB[i_std]/10)));
    }

    /* Bit Error Rate */
    float BER[num_EbN0dB];
    float BLER[num_EbN0dB];

/* Simulations */

    /*Simulations corresponding to each sig*/
    int i_sig = 0;
    int count[1];
    for (i_sig = 0; i_sig < num_EbN0dB; i_sig++)
    {
        printf("Count %d of %d\n", i_sig+1, num_EbN0dB);

        int err_count = 0;
        sig = sigma[i_sig];
        
        BER[i_sig] = 0;
        BLER[i_sig] = 0;

        /*No. of errors corresponding to each sig over num_sim simulations*/
        int num_err = 0;

        /*Each Simulations*/
        int i_num_sim = 0;
        for (i_num_sim = 0; i_num_sim < num_sim; i_num_sim++)
        {
            /*Generating a uniformly distributed binary random vector for message*/
            int msg[K];
            int i_msg = 0;
            for (i_msg = 0; i_msg < K - crc_l; i_msg++)
            {
                msg[i_msg] = 0*uni();
            }

            for (i_msg = K - crc_l; i_msg < K ; i_msg++)
            {
                msg[i_msg] = 0;
            }

            /* Encoding the Transmit vector */
            enc_start = clock();

            crcGen(msg, K, polynomial, crc_l);

            int u[N];
            int i_uN = 0;
            for (i_uN = 0; i_uN < N; i_uN++)
            {
                u[i_uN] = 0;
            }

            /*Assigning data to data indices*/
            int i_ud = 0;
            for (i_ud = 0; i_ud < K; i_ud++)
            {
                u[data_pos[i_ud]] = msg[i_ud];
            }

            encode(u, N);
            enc_end = clock();

            enc_cpu_time_used = enc_cpu_time_used + ( enc_end - enc_start );

            /*BPSK modulation and Adding Gaussian noise of zero mean and variance sigma^2*/
            /*Energy per bit*/
            int Eb = 1;
            int i_bpsk = 0;
            int x[N];
            float y[N];
            for (i_bpsk = 0; i_bpsk < N; i_bpsk++)
            {
                x[i_bpsk] = sqrt(Eb) * (1 - 2 * u[i_bpsk]);
                y[i_bpsk] = x[i_bpsk] + sig*randn(0,1);
            }

            /* Channel LLR calculation and Quantization of LLR values */
            float LLR[N];
            int LLR_Q[N];
            int i_ch = 0;
            int i_list = 0;
            for (i_ch = 0; i_ch < N; i_ch++)
            {
                LLR[i_ch] = 2*y[i_ch]/(pow(sig,2));
                // LLR[i_ch] = y[i_ch];
                // LLR_Q[i_ch] = floor((LLR[i_ch])/rmax*maxqr);
                // if ( LLR_Q[i_ch] > maxqr )
                // {
                //     LLR_Q[i_ch] = maxqr;
                // }
                // else if( LLR_Q[i_ch] < -(maxqr+1) )
                // {
                //     LLR_Q[i_ch] = -(maxqr+1);
                // }
            }

            /* Decoded message vector */
            int msg_cap[l*K];
            
            /* Successive Cancellation List Decoding */

            int crc_check[l];

            dec_start = clock();

            decode_unrolled(N, K, l, info_nodes, data_pos, LLR, msg_cap);

            crcDet(msg_cap, K, l, polynomial, crc_l, crc_check);

            int crc_ind = -1;
            int ii;
            for (ii = 0; ii < l; ii++)
            {
                if (crc_check[ii] == 0)
                {
                    crc_ind = ii;
                    break;
                }
            }

            if (crc_ind < 0)
            {
                crc_ind = 0;
            }

            dec_end = clock();

            dec_cpu_time_used = dec_cpu_time_used + ( dec_end - dec_start );

            int i_err = 0;
            int cur_err = 0;
            for (i_err = 0; i_err < K - crc_l; i_err++)
            {
                if (msg_cap[i_err + crc_ind*K] != msg[i_err])
                {
                    err_count++;
                    cur_err++;
                }
            }
            if (cur_err > 0)
            {
                BLER[i_sig]++;
            }
        } // end of no. of simulations loop

        printf("EbN0dB %0.2f\t", EbN0dB[i_sig]);
        BER[i_sig] = (double) err_count/(N*num_sim);
        BLER[i_sig] = (double) BLER[i_sig]/(num_sim);
        printf("BER %lf\t BLER %lf\n", BER[i_sig], BLER[i_sig]);

    } // end of noise var loop

	end = clock();
    cpu_time_used = ((float) (end - start)) / CLOCKS_PER_SEC;
    enc_cpu_time_used = (float) enc_cpu_time_used / CLOCKS_PER_SEC;
    dec_cpu_time_used_rec = ((float) dec_cpu_time_used_rec ) / CLOCKS_PER_SEC;
    dec_cpu_time_used = ((float) dec_cpu_time_used ) / CLOCKS_PER_SEC;
    printf("Time taken encode %d messages is %0.2f secs\n", num_sim*num_EbN0dB, enc_cpu_time_used);
    printf("Time taken decode %d codewords for unrolled decoder with l '%d' is %0.2f secs\n", num_sim*num_EbN0dB, l, dec_cpu_time_used);
    printf("Time taken to run %d simulations is %0.2f secs\n\n", num_sim, cpu_time_used);
    printf("Decoder throughput is %0.2f Mbps for list size %d\n", ((float)(num_sim*num_EbN0dB)/(dec_cpu_time_used*1000)),l);
    printf("Encoder throughput is %0.2f Mbps\n", ((float)(num_sim*num_EbN0dB)/(enc_cpu_time_used*1000)));

    int disp;
    printf("BER\n");
    for (disp = 0; disp < num_EbN0dB; disp++)
    {
        printf("%lf ", BER[disp]);
    }
    printf("\nBLER\n");
    for (disp = 0; disp < num_EbN0dB; disp++)
    {
        printf("%lf ", BLER[disp]);
    }

return 0;
}