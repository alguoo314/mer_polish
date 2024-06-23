#!/usr/bin/env python
import sys
import re
import os
import argparse
import math
import csv
from Bio.Align import PairwiseAligner
import dna_jellyfish as jf

def main(contigs,query_path,k,test,fix,fout,fixedout,db,thre,num_iter):
    try:
        divisor = 50
        qf  = jf.QueryMerFile(db)
        global solid_thre
        global debug 
        debug = False 
        global step
        step = max(2,round(k/8))
        global user_fix_choice
        user_fix_choice = fix
        solid_thre = thre #this is the threshold determined from the jellyfish histogram
        seq_dict = None
        for ite in range(num_iter+1): #num_iter rounds of fixing plus one more round to find the final q value
            query_path,seq_dict = iteration(num_iter,ite,qf,query_path,k,test,fix,fout,fixedout,db,seq_dict,divisor)
    except:
         exception_type, exception_object, exception_traceback = sys.exc_info()
         line_number = exception_traceback.tb_lineno
         print(line_number)
         print(sys.exc_info()) #to help debug                                  \                                                                   
         sys.exit(1)         
            

def iteration(num_iter,ite,qf,query_path,k,test,fix,fout,fixedout,database,seq_dict,divisor=50):   
    try:
        if ite == num_iter:
            fix=False
            fixedout = os.path.split(fixedout)
            fixedout = fixedout[0]+"_iter"+str(ite-1)+"_"+fixedout[1]
            
        if ite == 0:
            seq_dict = parse_fasta(query_path)
        wrong_kmers_count = 0
        fixed_bases_list = []
        total_wrong_kmers = 0
        total_kmers = 0
        fout = os.path.split(fout)
        fout = fout[0]+"_iter"+str(ite)+"_"+fout[1]
        for seqname,seq in seq_dict.items():
            total_kmers += len(seq)-k+1
            good_before = -1 #index of the last guaranteed good base before the mismatch
            i = 0 #first k mer at position 0
            wrong_kmers_count = 0                                                                                                       
            while i < len(seq)-k+1:
                mer_string = seq[i:k+i]
                N = mer_string.find('N') #ignore all kmers containing non acgt bases
                if N >= 0:
                    i+=(N+1)
                    continue
                n = mer_string.find('n')
                if n >= 0:
                    i+=(n+1)
                    continue
                match = re.match("^[ACTGactg]*$",mer_string) #other invalid characters other than N or n
                if match is None:
                    i +=1
                    continue
                    
                mer = jf.MerDNA(mer_string).get_canonical()
                occurrance = qf[mer]
                rolling_thre = 0
                if occurrance < solid_thre:
                    if debug:
                        print("Below normal threshold in " + seqname + " thresh " + str(solid_thre) + " position " +str(i))
                    i,seq,wrong_kmers_count,fixed_bases_list,break_while_loop = handle_bad_kmers(i,qf,seq,k,wrong_kmers_count,fix,fixed_bases_list,seqname,rolling_thre)
                    if break_while_loop:
                        break
                
                elif i > 0 and occurrance < qf[jf.MerDNA(seq[max(0,i-k):max(k,i)]).get_canonical()]/divisor: #above the solid threshold but > 1/divisor of the count of the kmer k bases before it
                    #check rolling average of previous k kmers
                    k_rolling_sum = 0
                    ind = max(0,i-k)
                    num=0
                    while ind < i:
                        num+=1
                        ind+=step
                        k_rolling_sum+=qf[jf.MerDNA(seq[ind:k+ind]).get_canonical()]
                    rolling_thre = round(k_rolling_sum/num/divisor)
                    if occurrance < rolling_thre:
                        if debug:
                            print("Below rolling threshold in " + seqname + " thresh " + str(rolling_thre) + " position " +str(i))
                        i,seq,wrong_kmers_count,fixed_bases_list,break_while_loop = handle_bad_kmers(i,qf,seq,k,wrong_kmers_count,fix,fixed_bases_list,seqname,round(k_rolling_sum/num/2))
                        if break_while_loop:
                            break
                    else:
                        i+=k-1
                    
                else: #good  kmer
                    i += k-1
            
                    
            seq_dict[seqname]=seq
            total_wrong_kmers += wrong_kmers_count


        if test == True:
            if ite == 0 or ite==num_iter:
                file_name = str(ite)+"qValCalcHelper.csv"
                with open(file_name,'a') as f:
                    f.write("{} {}\n".format(total_wrong_kmers,total_kmers))
                        
            
        if fix == True:
            base_fields = ['Contig', 'Base_coord', 'Original','Mutation']
            with open(fout,'w') as csvf:
                csvwriter = csv.writer(csvf,delimiter=' ')
                csvwriter.writerow(base_fields)
                csvwriter.writerows(fixed_bases_list)
        if ite == num_iter and user_fix_choice == True:
            i = 0
            with open(fixedout,'w') as of:
                for seqname,seq in seq_dict.items():
                    of.write(">{}\n".format(seqname))
                    new_seq = split_output(seq,60)
                    for l in new_seq:
                        of.write(l+"\n")
                    i+=1

        return fixedout,seq_dict
         
    except:
        exception_type, exception_object, exception_traceback = sys.exc_info()
        line_number = exception_traceback.tb_lineno
        print(line_number)
        print(sys.exc_info()) #to help debug
        sys.exit(1)




def split_output(seq, num_per_line=60): #make a new line after num_per_line bases
    lines = math.ceil(len(seq)/num_per_line)
    output = []
    for i in range(lines):
        output.append(seq[num_per_line*i:num_per_line*(i+1)])
    return output

         
def handle_bad_kmers(i,qf,seq,k,wrong_kmers_count,fix,fixed_bases_list,seqname,rolling_thre):
    thre = solid_thre
    if rolling_thre > 0:
       thre =  rolling_thre 
        
    j = i-1
    occurrance = qf[jf.MerDNA(seq[j:k+j]).get_canonical()] 
    while occurrance < thre and j>=0:
        j = j-1
        occurrance = qf[jf.MerDNA(seq[j:k+j]).get_canonical()]
    good_before = j+k-1 #the right base of a good kmer
    prev_good_count = qf[jf.MerDNA(seq[j:k+j]).get_canonical()]
    kmer_count = qf[jf.MerDNA(seq[i:k+i]).get_canonical()]                                                            

    if j == -1: #even the first kmer is bad                                                                                                                               
        good_before = -1

    if rolling_thre == 0: #solid threshold
        while kmer_count < thre and i < len(seq)-k+1: #i go forward
            i+=1
            kmer_count = qf[jf.MerDNA(seq[i:k+i]).get_canonical()]
    else:
        while kmer_count < thre and i < len(seq)-k+1: #rolling version only fix if number of bad kmers <= k
            if i-j > k:
                if debug:
                    print("the sequence of under-rolling-threshold kmers is longer than k")
                return i+1,seq,wrong_kmers_count,fixed_bases_list,False #ignore these kmers, move to the next index. In next iteration of the while loop, all k kmers before are the previously below-rolling-threshold kmers, which decreases the rolling-threshold
            i+=1
            kmer_count = qf[jf.MerDNA(seq[i:k+i]).get_canonical()]
    good_after = i #the first base of the first good kmer after the mismatch
    #A kmer is bad if it is below the threshold AND its count is less than 1/2 of the good k-mers before
    too_low_flag = False
    if (qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()] < solid_thre/2) and (qf[jf.MerDNA(seq[good_before-k+3:good_before+3]).get_canonical()] < solid_thre/2):
        too_low_flag = True #Below absolute thres = (relative) threshold/2, then it is an error kmer even if  it's > previous kemr count/2. Flag just for debugging purpose
    elif rolling_thre == 0: # this branch only apply to solid threshold scenerios
        while  qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()] >= prev_good_count/2 and good_before-k+1 < good_after: #kmer count > prev kmer count/2, then not a bad kmer.
            if good_before == -1:
                break
            if (prev_good_count >= thre/2 and (qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()] < thre/2) and (qf[jf.MerDNA(seq[good_before-k+3:good_before+3]).get_canonical()] < thre/2)):
                too_low_flag = True
                break
            #move to the next kmer
            prev_good_count = qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()]
            good_before +=1                        
        if good_before >= len(seq)-1:
            return i,seq,wrong_kmers_count,fixed_bases_list,True #break the while loop in the outer function. Switch to next seq
    #additional check for case 000...high high...000
    if max(0,good_before-k+2)+k+k >= len(seq):
        return max(0,good_before-k+2)+k+k,seq,wrong_kmers_count,fixed_bases_list,True #break the while loop in the outer function. Switch to next seq
    second = seq[max(0,good_before-k+2)+1:max(0,good_before-k+2)+k+1]
    k_minus_1 = seq[max(0,good_before-k+2)+k-2:max(0,good_before-k+2)+k+k-2]
    k_th =seq[max(0,good_before-k+2)+k-1:max(0,good_before-k+2)+k+k-1]
    k_plus_1 = seq[max(0,good_before-k+2)+k:max(0,good_before-k+2)+k+k]
    
    if qf[jf.MerDNA(second).get_canonical()] < thre and qf[jf.MerDNA(k_minus_1).get_canonical()] < thre and qf[jf.MerDNA(k_th).get_canonical()] < thre  and qf[jf.MerDNA(k_plus_1).get_canonical()]>=thre:
        good_after = max(0,good_before-k+2)+k #ie good_before+2
    to_be_fixed = seq[max(0,good_before-k+2):good_after+k-1]
    wrong_kmers_count += len([*range(max(0,good_before-k+2),good_after)])
    if debug:
        print("Found error in contig " + seqname + " " + to_be_fixed + " positions " + str(good_before) + " " +  str(good_after) + " bad kmers " + str(len([*range(max(0,good_before-k+2),good_after)])))
    if fix == True:
        if good_before < 0:
            return i,seq,wrong_kmers_count,fixed_bases_list,False 
        seq,fixed_base,original,fixed_ind = fixing_sid(seq,to_be_fixed,k,thre,qf,len([*range(max(0,good_before-k+2),good_after)]),good_before,good_after) #fix simple sub/insert/del cases
        if fixed_base != "nN":
            if rolling_thre > 0:
                if debug:
                    print("Rolling case FIXED!")
            if len(fixed_ind) == 1:
                fixed_bases_list.append([seqname,fixed_ind[0],fixed_base,original])
            else:
                fixed_bases_list.append([seqname,fixed_ind[0],fixed_base[0],original[0]])
                fixed_bases_list.append([seqname,fixed_ind[1],fixed_base[1],original[1]])
    return i,seq,wrong_kmers_count,fixed_bases_list,False 


def fixing_sid(seq,to_be_fixed,k,threshold,qf,num_below_thres_kmers,good_before,good_after): #fix sub and indel
    try:
        fixed_base = "nN"
        original = '-'
        fixed_ind = None

        if num_below_thres_kmers == k: #substitution or insertion
            b,fixed_subseq = fix_k_case_sub(to_be_fixed,k,threshold,qf)
            if b !=  None:
                original = "s"+seq[good_after-1]
                fixed_base = b
                fixed_ind = [good_after-1]
                seq = seq[:max(0,good_before-k+2)]+fixed_subseq+seq[good_after+k-1:]
            else:
                b,fixed_subseq = fix_insert(to_be_fixed,k,threshold,qf)
                if b != None:
                    original = "i"+seq[good_after-1]
                    fixed_base = '-'
                    fixed_ind = [good_after-1]
                    seq = seq[:max(0,good_before-k+2)]+fixed_subseq+seq[good_after+k-1:]       

        elif num_below_thres_kmers == k-1: #deletion or same base insertion or substitution
           removed_base,fixed_subseq = fix_del(to_be_fixed,k,threshold,qf)
           if removed_base != None: #deletion of a base
               original = "d-"
               fixed_ind = [good_after]
               seq = seq[:max(0,good_before-k+2)]+fixed_subseq+seq[good_after+k-1:]
               fixed_base = removed_base 
           else:
              left,right,l_or_r,fixed_subseq =  fixdiploid(to_be_fixed,k,threshold,qf,seq,good_before,good_after) #diploidy of two adjacent bases
              if l_or_r !=  None:
                  if l_or_r == "s": #lefting base is changed
                      original = "s"+seq[good_after-1]
                      fixed_base = str(left)
                      fixed_ind = [good_after-1]
                  else:#the righting base is changed
                      original = "s"+seq[good_before+1]
                      fixed_base = str(right)
                      fixed_ind = [good_before+1]
                  seq = seq[:max(0,good_before-k+2)]+fixed_subseq+seq[good_after+k-1:]
              else:
                  inserted_index,inserted_base,fixed_subseq = fix_same_base_insertion(to_be_fixed,k,threshold,qf,num_below_thres_kmers) #expensive, try last
                  if inserted_base != None:
                      original = "i"+inserted_base
                      fixed_base = "-"
                      seq = seq[:max(0,good_before-k+2)]+fixed_subseq+seq[good_after+k-1:]
                      fixed_ind = [inserted_index+max(0,good_before-k+2)] #inserted base at this index
 
        elif num_below_thres_kmers < k-1 and num_below_thres_kmers > 1 and len(to_be_fixed)>=k:#skip the good_before = -1 (ie first kmer is bad) case.
            removed_index,removed_base,fixed_subseq = fix_same_base_del(to_be_fixed,k,threshold,qf,num_below_thres_kmers) #the cheapest, try first
            if removed_base != None: #deletion of a base
                original = "d-"
                fixed_ind = [removed_index+max(0,good_before-k+2)]
                seq = seq[:max(0,good_before-k+2)]+fixed_subseq+seq[good_after+k-1:]
                fixed_base = removed_base
            else:
                left,right,l_or_r,fixed_subseq =  fixdiploid(to_be_fixed,k,threshold,qf,seq,good_before,good_after) #diploidy
                if l_or_r !=  None: 
                    if l_or_r == "s": #lefting base is changed
                        original = "s"+seq[good_after-1]
                        fixed_base = str(left)
                        fixed_ind = [good_after-1]
                    else:#the righting base is changed
                        original = "s"+seq[good_before+1]
                        fixed_base = str(right)
                        fixed_ind = [good_before+1]
                    seq = seq[:max(0,good_before-k+2)]+fixed_subseq+seq[good_after+k-1:]                
                else:
                    inserted_index,inserted_base,fixed_subseq = fix_same_base_insertion(to_be_fixed,k,threshold,qf,num_below_thres_kmers) #expensive, try last
                    if inserted_base != None:
                        original = "i"+inserted_base
                        seq = seq[:max(0,good_before-k+2)]+fixed_subseq+seq[good_after+k-1:]
                        fixed_base = "-"
                        fixed_ind = [inserted_index+max(0,good_before-k+2)] #inserted base is at this index

        elif num_below_thres_kmers > k: #two or more nearby errors.
            good_kmer_before = seq[good_before-k+1:good_before+1] 
            good_k_mer_after = seq[good_after:good_after+k] 
            fixed_seq = base_extension(len(to_be_fixed),qf,k,good_kmer_before,good_k_mer_after,threshold)
            if fixed_seq != None:
                fixed_ind = []
                fixed_base = []
                original = []
                aligner = PairwiseAligner()
                aligner.match_score = 0
                aligner.mismatch_score = -1
                aligner.open_gap_score = -1
                aligner.extend_gap_score = -1
                difference = aligner.align(fixed_seq, seq[good_before + 1:good_after])[0]
                fixed_seq_rep = difference[0]
                original_rep = difference[1]
                seq = seq[:good_before+1]+fixed_seq+seq[good_after:]
                for index in range(len(fixed_seq_rep)):
                    ori=original_rep[index]
                    changed = fixed_seq_rep[index]
                    if changed == ori:
                        continue
                    elif changed == "-": #fixed an insertion
                        fixed_base.append('-')
                        original.append("i"+ori)
                        fixed_ind.append(index+good_before+1)
                    elif ori == "-": #fixed a deletion
                        original.append("d-")
                        fixed_ind.append(index+(good_before+1))
                        fixed_base.append(changed)
                    else: #fixed a snp
                        original.append("s"+ori)
                        fixed_base.append(changed)
                        fixed_ind.append(index+(good_before+1))
                        
                        
        return seq,fixed_base, original, fixed_ind
    except:
        exception_type, exception_object, exception_traceback = sys.exc_info()
        line_number = exception_traceback.tb_lineno
        print(line_number)
        print(sys.exc_info()) #to help debug                                  \                                                                   
        sys.exit(1)         
                
def fixdiploid(seq_to_be_fixed,k,threshold,qf,full_seq,good_before,good_after):
    if debug:
        print("Trying diploid fix")
    #for fixing the first base of the last k-mer in L we also need to check the kmers before that contains this base.
    # for fixing the last b of first kmer in L we also need to check if the kmer starting at good_after is good still
    try:
        #if threshold > solid_thre: #We dont use complicated fixing method on kmers below the rolling threshold 
        #    return None,None,None,None
        left_bad = seq_to_be_fixed[len(seq_to_be_fixed)-k]
        right_bad = seq_to_be_fixed[k-1]
        left = left_bad
        right = right_bad
        good_before_starting_ind = max(0,good_before-k+1)
        base_after=''
        if good_after+k-1+int((k-1-len(seq_to_be_fixed)+k)/2) < len(full_seq):
            base_after = full_seq[good_after+k-1:good_after+k-1+int((k-1-len(seq_to_be_fixed)+k)/2)] #the end of the kmer starting with the right changed base = right index-left index -1 + good_after+k-1 
        else:
            base_after = full_seq[min(len(full_seq)-1,good_after+k-1):len(full_seq)]
        before_len=len(base_after)
        bases_before=full_seq[max(0,good_before_starting_ind-before_len+1):good_before_starting_ind+1]
        for x in 'ACTG':
            for y in 'ACTG':
                if x ==left_bad and y==right_bad:
                    continue
                if x!=left_bad and y!=right_bad:
                    continue
                trial = seq_to_be_fixed[:len(seq_to_be_fixed)-k]+x+seq_to_be_fixed[len(seq_to_be_fixed)-k+1:k-1]+y+seq_to_be_fixed[k:]
                check=bases_before+trial+base_after
                if check_sequence(check,qf,k,threshold):
                    left = x
                    right = y
                    if x == left_bad:
                        l_or_r = "e" 
                        
                    elif y == right_bad:
                        l_or_r = "s"
                        
                    else: #both changed not acceptable
                        continue
                    if debug:
                         print("Success diploid");
                    return(left,right,l_or_r,trial)
        return None,None,None,None
    except:
         exception_type, exception_object, exception_traceback = sys.exc_info()
         line_number = exception_traceback.tb_lineno
         print(line_number)
         print(sys.exc_info()) #to help debug                                  \                                                                   
         sys.exit(1)    
    


def fix_k_case_sub(seq_to_be_fixed,k,threshold,qf): #when number of conseuctive bad kmers is k
    bad_base = seq_to_be_fixed[k-1]
    if debug:
        print("Trying base substitution of " + bad_base + " in " +seq_to_be_fixed)
    for b in 'ACTG':
        trial = seq_to_be_fixed
        if b == bad_base:
            continue
        else:
            trial = trial[:k-1]+b+trial[k:]
            if check_sequence(trial,qf,k,threshold):
                if debug:
                    print("Success substitution "+ trial);
                return b,trial
    return None,None
                

def fix_insert(seq_to_be_fixed,k,threshold,qf):
    if debug:
        print("Trying to fix an insertion in " + seq_to_be_fixed)
    ind_to_be_removed = k-1
    base_to_be_removed = seq_to_be_fixed[ind_to_be_removed]
    seq_to_be_fixed = seq_to_be_fixed[:ind_to_be_removed] + seq_to_be_fixed[ind_to_be_removed+1:]
    if check_sequence(seq_to_be_fixed,qf,k,threshold):
            if debug:
                print("Success insertion "+seq_to_be_fixed)
            return base_to_be_removed,seq_to_be_fixed
    return None,None


def fix_del(seq_to_be_fixed,k,threshold,qf): 
    if debug:
        print("Trying to fix deletion in " + seq_to_be_fixed + " between " + seq_to_be_fixed[:k-1] + " and " +seq_to_be_fixed[k-1:])
    for alt in 'ATCG':
        trial = seq_to_be_fixed[:k-1]+alt+seq_to_be_fixed[k-1:]
        if check_sequence(trial,qf,k,threshold):
            if debug:
                print("Success deletion "+trial)
            return alt,trial
    return None,None


def fix_same_base_del(seq_to_be_fixed,k,threshold,qf,num_below_thres_kmers):
    if debug:
        print("Trying same base deletion in " + seq_to_be_fixed)
    if threshold > solid_thre: #We dont use complicated fixing method on kmers below the rolling threshold 
        return None,None,None
    sb = seq_to_be_fixed[k-2] #sb stands for same base
    fixed = False
    trial=seq_to_be_fixed
    new_bad = 0
    inserted = 0
    original_bad = len(seq_to_be_fixed)-k+1
    current_bad = original_bad
    max_insertions = original_bad
    if debug:
        print("Original bad "+str(original_bad) + " sb "+ sb)
    while inserted < max_insertions:
        new_bad=0
        trial = trial[:k-1]+sb+trial[k-1:]
        fixed = True
        inserted+=1
        for i in range(0,len(trial)-k+1):
            if qf[jf.MerDNA(trial[i:k+i]).get_canonical()] < threshold:
                fixed  = False
                new_bad +=1
        if fixed == True:
            if debug:
                print("Success same_base_deletion " + trial)
            return k-1,sb*inserted,trial
        if (new_bad >= current_bad):
            inserted = max_insertions
            break
        else: #added one base may have helped but need more
            current_bad = new_bad
            continue
    #now we try to insert a base before the first base of the first good k-mer
    if debug:
        print("Trying to insert a single base before the first good k-mer")
    for alt in 'ATCG':
        trial = seq_to_be_fixed[:k-2]+alt+seq_to_be_fixed[k-2:]
        if check_sequence(trial,qf,k,threshold):
            if debug:
                print("Success before_first_deletion " + trial)
            return k-2,alt,trial 
    return None,None,None

def fix_same_base_insertion(seq_to_be_fixed,k,threshold,qf,num_below_thres_kmers):
    if debug:
        print("Trying same base insertion in " + seq_to_be_fixed)
    if threshold > solid_thre: #We dont use complicated fixing method on kmers below the rolling threshold but above the solid (relative) threshold Rt
        return None,None,None
    ind_to_be_removed = k-1
    sb = seq_to_be_fixed[k-1] #sb stands for same base
    seq_to_be_fixed_local=seq_to_be_fixed
    fixed = False
    deleted = 0
    original_bad = len(seq_to_be_fixed)-k+1
    current_bad = original_bad
    max_deletions = original_bad
    if debug:
        print("Original bad "+str(original_bad) + " sb "+ sb + " num_kmers " +str(num_below_thres_kmers))
    while seq_to_be_fixed[k-1] == sb and deleted < max_deletions:
        current_bad -=1
        deleted +=1
        seq_to_be_fixed_local = seq_to_be_fixed_local[:ind_to_be_removed] + seq_to_be_fixed_local[ind_to_be_removed+1:]
        if len(seq_to_be_fixed_local) == k:
            break
        fixed=True
        new_bad = 0
        for i in  range(0,len(seq_to_be_fixed_local)-k+1):
            if qf[jf.MerDNA(seq_to_be_fixed_local[i:k+i]).get_canonical()] < threshold:
                fixed=False
                new_bad +=1
        if (fixed == True):
            if debug:
                print("Success same_base_insertion " + seq_to_be_fixed_local)
            return k-1,sb*deleted,seq_to_be_fixed_local
        if (new_bad >= current_bad):
            break
        else: #delete one more base
            current_bad = new_bad
            continue
    if debug:
        print("Let's try to delete a single base within a window of len(seq_to_be_fixed)-k and len(seq_to_be_fixed)")
    for i in range(len(seq_to_be_fixed)-k,len(seq_to_be_fixed)-1):
        trial = seq_to_be_fixed[:i]+seq_to_be_fixed[i+1:]
        #print("trying " +trial)
        if check_sequence(trial,qf,k,threshold):
            if debug:
                print("Success single_base_deletion " + trial + " i= " + str(i) + " Length = " + str(len(seq_to_be_fixed)))
            return i,seq_to_be_fixed[i],trial #Need further modification later because the index for the deleted base is different from case 1
    return None,None,None


def base_extension(len_seq_to_be_fixed,qf,k,good_kmer_before,good_k_mer_after,threshold):
    if len(good_kmer_before) < k or len(good_k_mer_after) < k or threshold > solid_thre:#We dont use complicated fixing method on kmers below the rolling threshold 
        return None
    bases = ["A", "C", "G", "T"]
    start_km1 = good_kmer_before[0:k-1]   # store the k-1 bases in a variable no need to carry these around
    min_overlap = 5
    for slack in range(10,11,1):
        paths = [] # array of all possible extensions
        max_ext = round((len_seq_to_be_fixed - 2*k)*1.2) + min_overlap + slack
        min_patch_len = round((len_seq_to_be_fixed - 2*k)/1.2) - slack                                                                                     
        paths.append(good_kmer_before[k-1:k])  # the last base of the initial k-mer makes the first path 
        if debug:
            print("Looking for a path "+str(len_seq_to_be_fixed) + " " + str(min_patch_len) +" "+ start_km1 + " " + paths[0] + " "+good_kmer_before+ " " +good_k_mer_after)

        for i in range(1,max_ext):
            paths = [l for l in paths if len(l) > 0]
            if len(paths) > 5000:
                if debug:
                    print("Too many paths")
                return None
            last_path = len(paths)  # since the number of paths will be changing we need to record this                                                                                       
            for p in range(last_path):
                if paths[p]=="":
                    continue
                km1 = (start_km1 + paths[p])[-k+1:]  # get the k-1 bases off the end                                                                         
                path_ext_count = 0  # this becomes 1 if we find an extension, and if no extension, then we delete the path
                for j in range(4):  # try to extend                                                       
                    score = qf[jf.MerDNA(km1 + bases[j]).get_canonical()]
                    if score >= threshold:
                        last_bases=km1 + bases[j]
                        if i >= min_overlap and i >= min_patch_len:
                            if last_bases[-min_overlap:]==good_k_mer_after[0:min_overlap]:
                                if path_ext_count: #check if this is the first extension
                                    path_connected=(start_km1 + (paths[p])[:-1] + bases[j] + good_k_mer_after[-(k-min_overlap):])[-(2*k-1):];
                                    return_path=((paths[p])[:-1] + bases[j])[1:-min_overlap]
                                else:
                                    path_connected=(start_km1 + paths[p] + bases[j] + good_k_mer_after[-(k-min_overlap):])[-(2*k-1):];
                                    return_path=(paths[p] + bases[j])[1:-min_overlap]
                                if debug:
                                    print("Candidate path "+ path_connected + " target " + good_k_mer_after + " iteration " +str(i))
                                if check_sequence(path_connected,qf,k,threshold):
                                    if i == min_overlap:
                                        if debug:
                                            print("Success path "+ start_km1 + paths[p]+ bases[j] + " slack " + str(slack) + " patch empty ")
                                        return None
                                    else:
                                        if debug:
                                            print("Success path "+ start_km1 + paths[p]+ bases[j] + " slack " + str(slack) + " patch "+ (paths[p]+bases[j])[1:-min_overlap])
                                        return return_path 
                        if path_ext_count == 0:  # first extension -- extend the current path                                                                                 
                             paths[p] += bases[j]
                             path_ext_count = 1
                        else:  # another extension -- add a new path
                             paths.append((paths[p])[:-1] + bases[j])
                if path_ext_count == 0:
                    paths[p] = ""  # no extensions, kill this path 
    return None

def check_sequence(trial,qf,k,threshold):
    fixed = True
    #we have to do this many times, so this function must be optimal
    #we first check the first and the last k-mer and then check in the middle skipping every "step" k-mers
    if qf[jf.MerDNA(trial[:k]).get_canonical()] < threshold:
        fixed = False
    else:
        if qf[jf.MerDNA(trial[-k:]).get_canonical()] < threshold:
            fixed = False
        else:
            for i in  range(step,len(trial)-k,step):
                if qf[jf.MerDNA(trial[i:k+i]).get_canonical()] < threshold:
                    fixed  = False
                    break
    return fixed

def count_bad_kmers_step(trial,qf,k,threshold):
    bad_kmers=0
    #we have to do this many times, so this function must be optimal
    #we first check the first and the last k-mer and then check in the middle skipping every "step" k-mers
    if qf[jf.MerDNA(trial[:k]).get_canonical()] < threshold:
        bad_kmers+=1
    if qf[jf.MerDNA(trial[-k:]).get_canonical()] < threshold:
        bad_kmers+=1
    for i in  range(step,len(trial)-k,step):
        if qf[jf.MerDNA(trial[i:k+i]).get_canonical()] < threshold:
            bad_kmers+=1
    return bad_kmers


def parse_fasta(query_file):
    f = open(query_file,"r")
    temp = ''
    seq = {}
    name = "placeholder"
    for line in f:
        if line.startswith(">"):
            seq[name] = temp
            name = line.split()[0][1:] #save the sequence name/number
            temp = ''
        else:
            temp += line.replace('\n','') #remove whitespaces in the sequence

    seq[name] = temp
    seq.pop("placeholder") #remove the first key put into the dict as a placeholder                                                                
    f.close()
    return seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default = None, help="The path to the .jf  database file. Not needed if --contigs is given.")
    parser.add_argument("--reads",nargs='+',default = None, help="The path to the .fasta file(s）containing the contigs to build the jellyfish database. Not needed if --db is provided")
    parser.add_argument("-q","--query", help = "The path to the .fasta query file")
    parser.add_argument("-thre","--threshold", type=int, default = None, help = "The threshold for an unreliable kmer.")
    parser.add_argument("-k","--ksize", type=int,help = "The kmer size")
    parser.add_argument("--test", action='store_true',help = "Ouput the total num of bad kmers, and an estimation for Q value")
    parser.add_argument("--fix", action='store_true', help="Output the index of fixed bases and output the new sequence")
    parser.add_argument("--fout",default = "fout.csv",help = "The path to output the index of the fixed bases." )
    parser.add_argument("-ff","--fixedfasta",default = "fixed_seq.fasta",help = "The path to output the fixed assembly sequences")
    parser.add_argument("-p","--passes", type=int, default = 2, help = "The number of iterations of fixing.")
    args = parser.parse_args()
    main(args.reads,args.query,args.ksize,args.test,args.fix,args.fout,args.fixedfasta,args.db,args.threshold,args.passes)
