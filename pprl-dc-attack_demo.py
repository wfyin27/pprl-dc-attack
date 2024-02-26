
"""
This is a demo version,a version capable of attacking 2sh will be released after the paper is published.
"""


import pandas as pd
from datasketch import MinHash, MinHashLSH
from time import *
import numpy as np
from collections import deque, defaultdict,Counter
import multiprocessing
import bitarray

epochs = 12

def loadCSVData(filePath):
    df = pd.read_csv(filePath,dtype=object,)# dtype=object指定每列都为字符型
    df = df.drop_duplicates(subset=None, keep='first', inplace=False)
    return df.values.tolist(),df.columns.values

# -----------------------------------------------------------------------------

def gen_dict_plain(rec_palin_list):
    """
    根据明文数据集构成明文的字典
    输入：
        rec_palin_list      明文数据列表
    
    输出：
        dit_same            没有变化的实体的明文记录 [[明文1],[明文2]]
        dit_part            发生改变的实体的明文记录 [[明文1],[明文2]]
        dit_none            单个实体的记录（一条记录）
    """
    dit = {}
    dit_change = {}
    dit_same = {}
    dit_none = {}
    for rec_lst in rec_palin_list:
        for r in rec_lst:   
            if r[0] not in dit:
                dit[r[0]] = [r[attr_index[0]:attr_index[-1]+1]]
            else:
                dit[r[0]].append(r[attr_index[0]:attr_index[-1]+1])

    for i in dit:
        if len(dit[i])==1:
            dit_none[i]=dit[i]
        elif dit[i][0]==dit[i][1]:
            dit_same[i]=dit[i]
        elif dit[i][0]!=dit[i][1]:
            dit_change[i]=dit[i]
    
    return dit_same,dit_change,dit_none   

# -----------------------------------------------------------------------------

def gen_dict_encode(rec_encode_list):
    """
    根据明文数据集构成明文的字典
    输入：
        rec_palin_list      密文数据列表
    
    输出：
        dit_same            没有变化的实体的密文记录 [[密文1],[密文2]]
        dit_part            发生改变的实体的密文记录 [[密文1],[密文2]]
        dit_none            单个实体的记录（一条记录）
    """
    dit = {}
    dit_change = {}
    dit_same = {}
    dit_none = {}
    for rec_lst in rec_encode_list:
        for r in rec_lst:  
            
            if r[0] not in dit:
                dit[r[0]] = [r[-1]]
            else:
                dit[r[0]].append(r[-1])
    for i in dit:
        if len(dit[i])==1:
            dit_none[i]=dit[i]
        elif dit[i][0]==dit[i][1]:
            dit_same[i]=dit[i]
        elif dit[i][0]!=dit[i][1]:
            dit_change[i]=dit[i]
            # s1 = bitarray.bitarray(dit_change[i][0])
            # s2 = bitarray.bitarray(dit_change[i][1])
            # print(list(s1.itersearch(1)))
            # print(list(s2.itersearch(1)))
            # s = s1 & s2
            # print(list(s.itersearch(1)))
            # if s1.count()-s.count()>s2.count()-s.count():
            #     dit_change[i][0],dit_change[i][1] = dit_change[i][1],dit_change[i][0]
            # exit()
    return dit_same,dit_change,dit_none

# -----------------------------------------------------------------------------

def get_plain_qgram(dit_plain_change):
    """
        提取变化记录的变化情况
        输入：
            dit_plain_change   所有变化的明文记录   [[明文1],[明文2]]
        输出：
            dit_q_gram_add     只增了qgram
            dit_q_gram_rm      只减少了qgram
            dit_q_gram_mix     既增加qgram又减少qgram

    """
    qm1=q-1
    dit_q_gram_add  = {}
    dit_q_gram_rm   = {}
    dit_q_gram_mix  = {}
    dits = {}
    for voter_id in dit_plain_change:

        attr_val1=dit_plain_change[voter_id][0]     #记录1的所有属性
        attr_val2=dit_plain_change[voter_id][1]     #记录2的所有属性
        # attr_val1  = ['james','james']
        # attr_val2 =  ['james','jarne']
        #对于RBF，应该考虑q-gram是在那个属性中，根据每个链接属性将属性值转为qgram集合
        attr_1_q_gram_set = []
        attr_2_q_gram_set = []

        if type(attr_val2[0])==float  :
            attr_val2[0]=''

        if type(attr_val2[1])==float  :
            attr_val2[1]=''

        if type(attr_val1[0])==float  :
            attr_val1[0]=''

        if type(attr_val1[1])==float  :
            attr_val1[1]=''

        for nums in range(len(attr_index)):
         
            attr_1_q_gram_set += [attr_val1[nums][i:i+q] + attr_symbol[nums] for i in range(len(attr_val1[nums]) - qm1)]
            attr_2_q_gram_set += [attr_val2[nums][i:i+q] + attr_symbol[nums] for i in range(len(attr_val2[nums]) - qm1)]
        attr_1_q_gram_set = set(attr_1_q_gram_set)
        attr_2_q_gram_set = set(attr_2_q_gram_set)
       
        # 根据qgram集合求变化前后的交集，即没有变化的qgram
        qgram_change    = attr_1_q_gram_set^attr_2_q_gram_set
   
      
        # 得到因为变化引起的，增加和减少的qgram
        if (attr_1_q_gram_set | qgram_change) == attr_1_q_gram_set:
            dit_q_gram_rm[voter_id]=list(qgram_change)
            # dit_q_gram_mix[voter_id+'$'] = [[len(qgram_change)]+list(qgram_change)+[voter_id+'$'],[0]+[]]
        elif (attr_2_q_gram_set | qgram_change) == attr_2_q_gram_set:
            dit_q_gram_add[voter_id] = list(qgram_change)
            # dit_q_gram_mix[voter_id+'$'] = [[len1]+[],[len(qgram_change)]+list(qgram_change)+[voter_id+'$']]
        else:
            rm      = list(qgram_change & attr_1_q_gram_set)
            add     = list(qgram_change & attr_2_q_gram_set)
    
            len1    = len(rm)
            len2    = len(add)
            if voter_id+'#' not in dit_bit_mix:
                continue
            dit_q_gram_mix[voter_id+'$'] = [[len1]+rm+[voter_id+'$'],[len2]+add+[voter_id+'$']]
            # print(voter_id+'$',len(rm),len(list(dit_bit_mix[voter_id+'#'][0].itersearch(1))))
            if len(rm) not in dits:
                dits[len(rm)] =[len(list(dit_bit_mix[voter_id+'#'][0].itersearch(1))),1]
            else:
                dits[len(rm)][0]+= len(list(dit_bit_mix[voter_id+'#'][0].itersearch(1)))
                dits[len(rm)][1]+= 1
                
    for i in dits:
        print(i,dits[i][0]/dits[i][1])
    return dit_q_gram_add,dit_q_gram_rm,dit_q_gram_mix

# -----------------------------------------------------------------------------

def get_enc_bit(dit_enc_cahnge):
    """
        提取变化记录的变化情况
        输入：
            dit_enc_cahnge   所有变化的明文记录   [[密文1],[密文2]]
        输出：
            dit_bit_add     bit 0->1
            dit_bit_rm      bit 1->0
            dit_bit_mix     混合 

    """
    
    dit_bit_add = {}
    dit_bit_rm  = {}
    dit_bit_mix = {}

    for i in dit_enc_cahnge:

        bitarray_a = bitarray.bitarray(dit_enc_cahnge[i][0])
        bitarray_b = bitarray.bitarray(dit_enc_cahnge[i][1])
    
        s = bitarray_a & bitarray_b

        # 使用 ^ 运算符计算差集
        difference_rm = s ^ bitarray_a
        difference_add = s ^ bitarray_b
        # if difference_rm.count()<difference_add.count():
        #     difference_rm,difference_add=difference_add,difference_rm
        
        if difference_add.count()==0 and difference_rm.count()>0:  
            dit_bit_rm[i]=difference_rm
        elif difference_rm.count()==0 and difference_add.count()>0:  
            dit_bit_add[i]=difference_add
        else:
            dit_bit_mix[i+'#']=[difference_rm,difference_add,[i+'#'],0,0,0,0,0,0]
        
    
    return dit_bit_add,dit_bit_rm,dit_bit_mix


def minlsh(dit,dit2):
    block ={}
    lsh = MinHashLSH(threshold=0.4, num_perm=64)
    m1 = MinHash(num_perm=64)
   
    for i in dit2:
   
        set1 =set(dit2[i])
        for d in set1:
            m1.update(d.encode('utf8'))

        # Create LSH index
        lsh.insert(i, m1)
        m1.clear()

    m2 = MinHash(num_perm=64)
    for i in dit:
        set1 =set(dit[i][0])
        for d in set1:
            m2.update(d.encode('utf8'))
        result = lsh.query(m2)
        result=list(set(result))
        block[i]=result
        m2.clear()
    return block

def minlsh_qgram(dit_qgram_mix):
    block ={}
    block2= {}
    block3 = {}
    lsh = MinHashLSH(threshold=0.02, num_perm=64)
    m1 = MinHash(num_perm=64)

    lsh2 = MinHashLSH(threshold=0.02, num_perm=64)
    for i in dit_qgram_mix:
        rm  = dit_qgram_mix[i][0][1:-1]
        add = dit_qgram_mix[i][1][1:-1]  
        set1 =set(rm)
        for d in set1:
            m1.update(d.encode('utf8'))
        # Create LSH index
        lsh.insert(i, m1)
        m1.clear()

        set2 = set(add)
        for d in set2:
            m1.update(d.encode('utf8'))
        # Create LSH index
        lsh2.insert(i, m1)
        m1.clear()

    m3 = MinHash(num_perm=64)
    for i in dit_qgram_mix:
        rm  = dit_qgram_mix[i][0][1:-1]
        add = dit_qgram_mix[i][1][1:-1] 
        set1 =set(rm)
        for d in set1:
            m3.update(d.encode('utf8'))
        result = lsh.query(m3)
        result=list(result)
        block[i]=result

        result2 = lsh2.query(m3)
        result2=list(result2)
        block2[i]=result
        m3.clear()


        set1 =set(add)
        for d in set1:
            m3.update(d.encode('utf8'))
        result3 = lsh2.query(m3)
        result3=list(result3)
        block3[i]=result
        m3.clear()
    return block,block2,block3
import cosinlsh

def compare_minlsh(dit,dit2):
    """
    Optionally, use locally sensitive hash clustering or feature value based clustering
    """

    #-----------------------------use locally sensitive hash clustering-------------------------

    # lenght = 0
    # for i in dit2:
    #     lenght=max(lenght,dit2[i][1]+dit2[i][3]+dit2[i][5])
    # dit_qgram = {}
    # dit_bit = {}
    # for i in dit:

    #     result =dit[i][0]+ dit[i][2] + dit[i][4] # 取前 m 个元素
    #     result =result[:lenght]
    #     result += [0] * (lenght - len(result))  # 如果不足 m 个元素，用 0 补齐
    #     dit_bit[i]=result

    # for i in dit2:
    #     result =dit2[i][0]+ dit2[i][2] + dit2[i][4] # 取前 m 个元素
    #     result =result[:lenght]
    #     result += [0] * (lenght - len(result))  # 如果不足 m 个元素，用 0 补齐
    #     dit_qgram[i]=result 


    # Graph_sim_hash = cosinlsh.CosineLSH(lenght, 500,17)
    # qg_sim_hash_dict = Graph_sim_hash.gen_sim_hash(dit_qgram)
    # ba_sim_hash_dict = Graph_sim_hash.gen_sim_hash(dit_bit)

    # plain_sim_block_dict, encode_sim_block_dict = \
    #     Graph_sim_hash.hlsh_blocking(qg_sim_hash_dict,
    #                                     ba_sim_hash_dict,
    #                                     10,
    #                                     2,
    #                                     500,
    #                                     17)
    # dit_block = {}
    # for p in plain_sim_block_dict:
    #     for m in plain_sim_block_dict[p]:
    #         if m not in dit_block:
    #             dit_block[m]=[p]
    #         else:
    #             dit_block[m].append(p)
    #         print(len(dit_block[m]),len(plain_sim_block_dict))
    # block = {}
    # k=0
    # n=0
    # for p in dit_block:
    #     s=[]
    #     for m in dit_block[p]:
    #         s+=list(encode_sim_block_dict[m])
        
    #     block[p]=set(s)
    #     if p[:-1]+'#' in block[p]:
    #         k+=1
    #     n+=len(block[p])
    # print(k,len(dit_block),n//len(dit_block))

    # exit()
    # print(encode_sim_block_dict)

    #-----------------------------feature value based clustering-------------------------

    block ={}
    lsh = MinHashLSH(threshold=0.35, num_perm=64)
    m1 = MinHash(num_perm=64)
   
    for i in dit:
        s1 = []
        for j in dit[i][0]:
            s1.append(str(round(j / 20) * 20)+'-0')
      
        for j in dit[i][2]:
            s1.append(str(round(j / 20) * 20)+'-2')
          
        for j in dit[i][4]: 
            s1.append(str(round(j / 20) * 20)+'-4')
          
        # s1.append(s)
 
      
        set1 =set(s1)
        for d in set1:
            m1.update(str(d).encode('utf8'))

        # Create LSH index
        lsh.insert(i, m1)
        m1.clear()

    m2 = MinHash(num_perm=64)
    k=0
    for i in dit2:
        s1 = []
        for j in dit2[i][0]:
            s1.append(str(round(j / 20) * 20)+'-0')
          
        for j in dit2[i][2]:
            s1.append(str(round(j / 20) * 20)+'-2')
         
        for j in dit2[i][4]:
            s1.append(str(round(j / 20) * 20)+'-4')
           
        set1 =set(s1)
        for d in set1:
            m2.update(str(d).encode('utf8'))
        result = lsh.query(m2)
        result=list(set(result))
        if i[:-1]+'#' in result:
            k+=1
        block[i]=result
        m2.clear()
    print(k)
    return block





def get_feature_rm(dit):
    begin_time = time()

    block = {}
    for i in dit:
        dit_feature = {}
        for k in dit:
            s = list(set(dit[i]).intersection(dit[k]))
            l = len(s)
            if l > 0  and i!=k:
                if l not in block:
                    dit_feature[l]=[s]
                else:
                    dit_feature[l].append(s)
        block[i]=dit_feature
    end_time     =   time()
    print('end_time',end_time-begin_time)
   


#=========================构建q-gram图
def gen_qgram_grap(dit_qgram_mix):
    """
        构建qgram产生每个实体的qgram变化特征

        输入:
            dit_qgram_mix       每个实体qgram的变化情况  [[减少的qgram],[增加的qgram]]
        
        输出：
            dit_qgram_mix         更新dit_qgram_mix,添加每个实体变化的特征
            len_qgram_rm          减少的qgram中,具有相同数量变化的分布
            len_qgram_add         增加的qgram中,具有相同数量变化的分布
            len_qgram_mix         减少-增加的qgram中,具有相同数量变化的分布
    """

    #------------------step1 统计每个q_gram对应有那些实体 ，dit[q_gram]=[实体1，实体2]  即那些实体都改变了该qgram--------------------

    total_dit_rm    = defaultdict(lambda:[])          # 
    total_dit_add   = defaultdict(lambda:[])

    freq_one_qgram_rm  = defaultdict(lambda:0)
    freq_one_qgram_add  = defaultdict(lambda:0)    

    for i in dit_qgram_mix:
        rm  = dit_qgram_mix[i][0]
        add = dit_qgram_mix[i][1] 
        for j in rm[1:-1]:
            freq_one_qgram_rm[j]+=1
            total_dit_rm[j].append(rm[-1])
            
        for j in add[1:-1]:
            freq_one_qgram_add[j]+=1
            total_dit_add[j].append(add[-1])

    #-------------------------step2 构建q_gram图--------------------------------------------

    paramater = [   [total_dit_rm,freq_one_qgram_rm,'rm'],\
                    [total_dit_add,freq_one_qgram_add,'add'],\
                    [total_dit_add,freq_one_qgram_add,'mix']]
    
    keys_list = list(dit_qgram_mix.keys())
    epoch = epochs
    batch_size = int(len(dit_qgram_mix)/epoch)

    for i in paramater:
        arg1 = dict(i[0])
        arg2 = dict(i[1])
        arg4 = i[2]
        # gen_qgram_frep_grap(dit_qgram_mix,dit_qgram_mix,arg1,arg2,arg3,arg4)
        with multiprocessing.Pool(3) as  pool:
            for batch in batcher_iter(keys_list, batch_size):
             
                this_result =  pool.apply_async(gen_qgram_frep_grap,(batch,dit_qgram_mix,arg1,arg2,arg4),callback=save_result_ferp_grap)

            pool.close()
            pool.join()
        fit = defaultdict(lambda:0)
        for j in len_qgram_change[arg4]:
            for k in j:
                fit[k]+=j[k]
        len_qgram_change[arg4] = dict(fit)

def save_result_ferp_grap(result):

    len_qgram_change[result[2]]+=[result[0]]

    for i in result[1]:
        dit_qgram_mix[i].append( result[1][i])
 
 
def gen_qgram_frep_grap(batch,dit_qgram_mix,total_dit_change,freq_one_qgram,option):

    dit = {}

    len_qgram_change    = defaultdict(lambda:0)
    for i in batch:
        # 什么变化的qgram
        if option=='rm' or option=='mix':
            change_qgram  = dit_qgram_mix[i][0]
        elif option=='add':
            change_qgram   = dit_qgram_mix[i][1] 

        #-------------------------step2-1 统计每个实体和其他实体具有类似qgram变化的情况 --------------------------------------------

        # 统计和其他实体具有相同q_gram变化的数量
        s_change = []
        lst_change = defaultdict(lambda:[])
        
        for j in change_qgram[1:-1]:
            if j in total_dit_change:
                s_change += total_dit_change[j]
        s_change    = Counter(s_change)

        for j in change_qgram[1:-1]:
            if j in total_dit_change:
                for k in total_dit_change[j]:
                    if s_change[k]==1:
                        continue
                    lst_change[k].append(j)
            
        # 删掉自己和自己
        if len(change_qgram[1:-1])>1 and option!='mix': 
            del lst_change[i]

        # 只计算多个q_gram同时相似，加快速度
        lst  = gen_feature_qgram(lst_change,len_qgram_change)
        
        # 补上单个q_gram相同的变化
        for i2 in change_qgram[1:-1]:
            if i2 not in freq_one_qgram:
                continue
            lst[i2]  = freq_one_qgram[i2] 
            if option!='mix':
                 lst[i2]-=1
            for j in lst:
                if i2 in j and i2!=j:
                    lst[i2]-=lst[j]
            len_qgram_change[1]+=lst[i2]
        # print(option,lst)
        # exit()
   
        dit[i]=lst

    len_qgram_change = dict(sorted(len_qgram_change.items(),key=lambda x:x[0],reverse=False))
    # print('asdasdsad')
    # print(f'qgram{option}的总数：{len_qgram_change}')

    return [len_qgram_change,dit,option]   #用于构建转换表



def gen_feature_qgram(lst,len_qgram):
    """
        统计每个实体qgram变化时,和其他实体具有相同qgram变化的频率
        输入：
            lst:   一个字典 表示哪些实体和目标实体变化了某些相同的q_gram
            rm:     目标实体的qgram变化情况
        输出：
            len_qgram:  一个列表 统计每个实体与其他实体相同qgram变化的总数
            lst_rm:     一个字典  统计了和目标实体相同qgram变化的频率   例如dit[entity_id] = {'qgram1#qgram2#':22} 有22个实体都同时变化了qgram1和qgram2
    """
    lst_rm = defaultdict(lambda:0)
 
    for i in lst:
        feature =lst[i]
        
        feature.sort()
        len_qgram[len(feature)] += 1

        feature = ''.join(feature)
        lst_rm[feature] += 1

    return dict(lst_rm)

#=========================构建bit图           
def gen_bit_grap(dit_bit_mix,len_qgram_rm,len_qgram_add,len_qgram_mix):
#0 for i in range(len(dit_bit_mix))
    dit_rm = defaultdict(lambda:np.array([0 for i in range(len(dit_bit_mix))], dtype=np.uint8))
    dit_add= defaultdict(lambda:np.array([0 for i in range(len(dit_bit_mix))], dtype=np.uint8))
    dit_rm_bitarray = {}
    dit_add_bitarray= {}   
    dit_rel_id_nums = []     # 记录id和nums的映射

    k=0
    for i in dit_bit_mix:
        rm  = dit_bit_mix[i][0]
        add = dit_bit_mix[i][1]
        for j in list(rm.itersearch(1)):
            dit_rm[j][k]=1

        for j in list(add.itersearch(1)):
            dit_add[j][k]=1
        dit_rel_id_nums.append(i)
        k+=1
        dit_rm_bitarray[i]  =   rm
        dit_add_bitarray[i] =   add

    dit_rm=dict(dit_rm)
    dit_add=dict(dit_add)

    begin_time =   time()
    #-------------------------需要丢弃那些长度的变化，以及估计大于2个变化的长度是多少

    key_drop_rm,more_2_rm,table_bit_rm      = gen_key_drop_two(dit_rm_bitarray,dit_rm_bitarray,'rm',len_qgram_rm,dit_rm)
    key_drop_add,more_2_add,table_bit_add   = gen_key_drop_two(dit_add_bitarray,dit_add_bitarray,'add',len_qgram_add,dit_add)
    key_drop_mix,more_2_mix,table_bit_mix   = gen_key_drop_two(dit_rm_bitarray,dit_add_bitarray,'mix',len_qgram_mix,dit_add)

    key_drop_rm=max(key_drop_rm,key_drop_add,key_drop_mix)
    end_time =   time()
    print('产生全局key的时间:',end_time-begin_time)
 
    #--------------------------------生成每种变化的，每个实体的变化特征----------------------------------
    begin_time =   time()
    get_feature(dit_rm_bitarray,dit_rm_bitarray,'rm',0,global_a,table_bit_rm,key_drop_rm,more_2_rm,dit_rel_id_nums)
    get_feature(dit_add_bitarray,dit_add_bitarray,'add',1,global_a,table_bit_add,key_drop_rm,more_2_rm,dit_rel_id_nums)
    get_feature(dit_rm_bitarray,dit_add_bitarray,'mix',2,global_a,table_bit_mix,key_drop_rm,more_2_mix,dit_rel_id_nums)

    end_time =   time()
    print('julei:',end_time-begin_time)
  
#---------------------------------------------------------------------

def gen_key_drop(batch,dit_bitarray1,dit_bitarray2,option,dit_rm):
   
    len_bit2 = defaultdict(lambda:0)
    tmp_rel_id = {}
    for i in  batch:
        tmp = []
        rm = dit_bitarray1[i]
        for j in list(rm.itersearch(1)):
            if j not in dit_rm:
                continue
            if len(tmp)==0:
                tmp=np.copy(dit_rm[j])
            else:
                tmp+=dit_rm[j]
        tmp_rel_id[i] = tmp
        tmp = tmp[tmp != 0]
 
        unique_numbers, counts = np.unique(tmp, return_counts=True)
        for num, count in zip(unique_numbers, counts):
            len_bit2[num]+=count
     
        if option!='mix':
            len_bit2[rm.count()]-=1
  

    len_bit= {}
    for i in len_bit2:
        if len_bit2[i]!=0 and i !=0:
            len_bit[i]=len_bit2[i]
   
    return [option,len_bit,tmp_rel_id]



#---------------------------------------------------------------------
def gen_key_drop_two(dit_bitarray1,dit_bitarray2,option,len_qgram,dit_rm):

    keys_list = list(dit_bitarray1.keys())
    epoch = epochs
    batch_size = int(len(dit_bitarray1)/epoch)

    #gen_key_drop(keys_list,dit_bitarray1,dit_bitarray2,option,dit_rm)
    with multiprocessing.Pool(epoch) as  pool:
        for batch in batcher_iter(keys_list, batch_size):
      
            this_result =  pool.apply_async(gen_key_drop,(batch,dit_bitarray1,dit_bitarray2,option,dit_rm),callback=save_result_drop_key)

        pool.close()
        pool.join()

    len_bit = {}
    for i in drop_key_len[option]:
        for j in i:
            if j not in len_bit:
                len_bit[j] = i[j]
            else:
                len_bit[j] += i[j]

    key_drop,more_2=  find_nums_bit_drop(len_bit,len_qgram)
    print(f"{option}丢弃长度为{key_drop}的相似bit变化,大于2个变化的位{more_2}")

    for i in len_bit:
         if i<=key_drop:
            len_bit[i]=0
    table_bit = gen_transform_table(len_qgram,len_bit)
    return key_drop,more_2,table_bit

def save_result_drop_key(result):
    drop_key_len[result[0]].append(result[1])
    for i in result[2]:
       
        rel_id[result[0]][i]=result[2][i]


def save_result(result):
    dit_enc_feature[result[1]].update(result[0])

def tmp2(batch,dit_bitarray1,dit_bitarray2,global_a,key_drop,more_2,keys,table_bit,option,dit_rel_id_nums):
    tmp_lst_feature = {}

    for i in  batch:
        rm = dit_bitarray1[i]
        lst = []
        lst2 = []


        indices = np.where(rel_id[option][i] > key_drop)[0]

        for j in indices:
            if dit_rel_id_nums[j]==i and option!='mix':
                continue
            rm2 = dit_bitarray2[dit_rel_id_nums[j]]
            s = rm&rm2
            lst.append(s)
            lst2.append(s.count())
   

        

        lst_rm,tmp_rm = gen_feature_bit(lst,key_drop,more_2,global_a,len(dit_bitarray2),lst2) 
        f_rm_1,f_rm_2,max_bit_rm,a,tmp_bit_dit = gen_feature(lst_rm,table_bit,tmp_rm,global_a)
        tmp_lst_feature[i]=[f_rm_1,f_rm_2,max_bit_rm,a]
       
  
    return tmp_lst_feature,keys

def get_feature(dit_bitarray1,dit_bitarray2,option,keys,global_a,table_bit,key_drop,more_2,dit_rel_id_nums):

    keys_list = list(dit_bitarray1.keys())
    epoch = epochs
    batch_size = int(len(dit_bitarray1)/epoch)

    #tmp2(dit_bitarray1,dit_bitarray1,dit_bitarray2,global_a,key_drop,more_2,keys,table_bit,option,dit_rel_id_nums)
    with multiprocessing.Pool(epoch) as  pool:
        for batch in batcher_iter(keys_list, batch_size):
            this_result =  pool.apply_async(tmp2,(batch,dit_bitarray1,dit_bitarray2,global_a,key_drop,more_2,keys,table_bit,option,dit_rel_id_nums),callback=save_result)

        pool.close()
        pool.join()

    

def find_nums_bit_drop(len_bit,len_qgram):
    """
    由于存在哈希冲突,需要判断一个qgram变化影响至少应该有多少个bit位表示
    输入：
    len_bit : bit位变化的频率分布
    len_qgram : qgram变化的频率分布

    输出：
        至少几位bit位
    """

    len_bit = dict(sorted(len_bit.items(),key=lambda x:x[0],reverse=False))
   
    sun_qgram_nums = sum(len_qgram.values())
    sun_qgram_1 = sun_qgram_nums-len_qgram[1]
    sun_bit_nums = sum(len_bit.values())
    tmps = 0
    nums_drop_min_nums = 0
    # print('asdas',sun_qgram_nums)
    # print(len_bit)
    # print(len_qgram)
    # print('asd',sun_bit_nums)
    nums_drop_min = 0
    more_2 = 0
    for i in len_bit:
        tmps+=len_bit[i]
     
        if sun_bit_nums-tmps<sun_qgram_1*0.8:           
            more_2= nums_drop_min
            break
        nums_drop_min =i
    sun_qgram_nums=sun_qgram_nums*1.25
    tmps = 0
    for i in len_bit:
        tmps+=len_bit[i]
        nums_drop_max_nums = sun_bit_nums-tmps
        # print(i,tmps,sun_bit_nums,nums_drop_max_nums,sun_qgram_nums)
        if nums_drop_max_nums<sun_qgram_nums and nums_drop_min_nums>=sun_qgram_nums:
         
            if abs(nums_drop_min_nums-sun_qgram_nums)>=abs(nums_drop_max_nums-sun_qgram_nums):
                return i,more_2
            else:
                return nums_drop_min,more_2
        nums_drop_min =i
        nums_drop_min_nums = sun_bit_nums-tmps

    return 0,more_2


cc = 0
def gen_feature_bit(lst_same_bit,key_drop,more_2,global_a,len_dit_bitarray,lst2): #more_2,option,i
    tmp = {}
    lst =  {}
    block = 0
    data =[]
    tmpbit = bitarray.bitarray(1000)    
    tmpbit.setall(1)
    tmp[block] = [0,tmpbit,0,0]
    lst[block]=[0,0]

    dit_block = defaultdict(lambda:0)
    sorted_indices =  sorted(range(len(lst2)), key=lambda k:lst2[k],reverse=True)
    lst2.sort(reverse=True)

    # t=0
    left=0



    for index, key in enumerate(lst2):
        if key<=more_2:
            left=index
            break
    for j in range(left,len(sorted_indices)):
        same_atr =lst_same_bit[sorted_indices[j]]
        indexes = list(same_atr.itersearch(1))
        # print(indexes)
        max_cluster = 0
        max_index = 0
        
        # tmplsy = [0 for i in range(block+1)]
    
        # for i in range(0,min(len(indexes),5)):
        #     tmplsy[dit_block[indexes[i]]]+=1
        # tmplsy =  sorted(range(len(tmplsy)), key=lambda k:tmplsy[k],reverse=True)
        # tmplsy=tmplsy[:min(3,block)]
        # tmplsy = []
        # # max_index= tmplsy.index(max(tmplsy))
        # for i in range(0,min(len(indexes),5)):
        #     tmplsy.append(dit_block[indexes[i]])
    
        # tmplsy = set(tmplsy)
        # max_cluster=(tmp[max_index][1] & same_atr).count()/tmp[max_index][1].count()
        
        # if tmp_same_atr!=0 and (tmp_same_atr & same_atr).count()/tmp_same_atr.count()>0.8:
        #     max_index = tmp_index
        #     max_cluster=(tmp[max_index][1] & same_atr).count()/tmp[max_index][1].count()
        
        # else:
        
        for k in lst:
            s=(tmp[k][1] & same_atr).count()/tmp[k][1].count()
            
            if max_cluster<s :
                max_cluster = s
                max_index = k
            if s>0.8:
                break
        
        if max_cluster >=0.55:
            lst[max_index][0]+=1
            lst[max_index][1]=max(lst[max_index][1],same_atr.count())
            tmp[max_index][2][sorted_indices[j]]=1 
        
            # if lst[max_index][0]>300:
            #     continue
            for k3 in indexes:
                tmp[max_index][0][k3]+=1
            if lst[max_index][0]%5!=0:
                continue
            tmp[max_index][1].setall(0)
            for k4 in  tmp[max_index][0]:                              
                if tmp[max_index][0][k4]/lst[max_index][0]>=thor_signature:
                    tmp[max_index][1][k4]=True
                    dit_block[k4]=max_index
        else:
            # if len(tmp)>max_f_rm:
            #     continue

            lst[block]=[1,len(indexes)]
            tmp_origin =defaultdict(lambda:0)
            for t2 in indexes:
                tmp_origin[t2] = 1
                dit_block[t2]=block
            tmpbit = bitarray.bitarray(len_dit_bitarray)    
            tmpbit.setall(0)
            tmpbit[sorted_indices[j]]=1
            tmp[block] =[tmp_origin,same_atr,tmpbit]
            block+=1

    keys_list = list(tmp.keys())
    for i2 in keys_list:
     
        if tmp[i2][1].count()<=key_drop  or lst[i2][0]<global_a:
            del lst[i2]
            del tmp[i2]
    
        
    for i2 in lst:
        for j in lst:
            if i2==j or lst[j][0]==0 or lst[i2][0]==0:
                continue
            if (tmp[i2][1]&tmp[j][1]).count()==tmp[j][1].count() :
                lst[j][0]+=lst[i2][0]
                lst[i2][0]=0  
                tmp[j][2]|=tmp[i2][2]
                break
    # print(len(lst),len(lst_same_bit))
    for i2 in lst:
        if lst[i2][0]==0:
            continue 
        for k in range(len(lst_same_bit)):
            j=sorted_indices[k]
            if tmp[i2][2][j]==False:
                if (tmp[i2][1]&lst_same_bit[j]).count()/(tmp[i2][1]).count()>=thor_signature:
                    lst[i2][0]+=1
                    lst[i2][1] = max(lst[i2][1],lst_same_bit[j].count())
            
    # if i[0:-1] =='AK61189':
    #     print(data)
    #     exit()
    # print(t,len(lst_same_bit))
        
    return lst,tmp



#转换q-gram变化个数和bit为变化个数
def gen_transform_table(len_qgram,len_bit):
    len_qgram = dict(sorted(len_qgram.items(),key=lambda x:x[0]))
 
    len_bit = dict(sorted(len_bit.items(),key=lambda x:x[0]))
    len_qgram = count_dit(len_qgram)
    len_bit = count_dit(len_bit)
  
    table_bit_rm = {}
    tmp = 1
    for i in len_bit:
        for k in len_qgram:
            if len_bit[i]<=len_qgram[k]:
                table_bit_rm[i]=[k]
                if len_qgram[k]-len_qgram[tmp]<0.01 and tmp !=k:
                    table_bit_rm[i].append(tmp)
                break
            tmp = k
    gap =3
    tmp_table_bit_rm = {}
    for i in table_bit_rm:
        tmp_table_bit_rm[i] =table_bit_rm[i] 
        for z in range(i-gap,i+gap):
            if z in table_bit_rm:
                tmp_table_bit_rm[i]=list(set(tmp_table_bit_rm[i]).union(set(table_bit_rm[z])))
               
    return tmp_table_bit_rm

def count_dit(di1):
    tmp = 0
    m =0
    for i in di1:
        di1[i]+=tmp
        tmp= di1[i]
        m = di1[i]
    for i in di1:
        di1[i]/=m
    return di1

def count_elements_greater_than(arr, target):
    """
    计算数字数组中大于目标数字的元素个数
    参数：
    - arr: 数字数组
    - target: 目标数字
    返回值：
    - 符合条件的元素个数
    """
    # 将数字数组转换为 NumPy 数组
    arr_np = np.array(arr)
    # 使用 NumPy 的条件语句生成一个布尔数组，表示数组中每个元素是否大于目标数字
    bool_array = arr_np > target
    # 使用布尔数组进行索引，得到符合条件的元素，并计算个数
    count = np.sum(bool_array)
    # 返回符合条件的元素个数
    return count


def get_gobal_key(dit_qgram_mix):
    """
    对明文特则统计分析,产生一个全局key,小于这个key值频率特征将被丢弃,允许平均每个实体最多丢弃1个特征
    输入 :
        dit_qgram_mix: 明文特征
    输出：
        key1 : 计算频率为key1以上的特征值
        key2 : 计算相似度时，给每个属性的权重,为删去不明显特征后，平均每个实体的特征数量
    """
    dits_rm = []
    dits_add = []
    dits_mix = []
    for i in dit_qgram_mix:
        # print(dit_qgram_mix[i])
        # exit()
        dits1 = get_qgram_nums(dit_qgram_mix[i][2])
        dits2 = get_qgram_nums(dit_qgram_mix[i][3])
        dits3 = get_qgram_nums(dit_qgram_mix[i][4])
        # if len(dits1)==0:
        #     continue
        dits_rm+=dits1
        dits_add+=dits2
        dits_mix+=dits3
    # print(min(dits_rm),sum(dits_rm)/len(dit_qgram_mix))
    # exit()
    # for i in dits_rm:
    #     print(dits_rm[i])
    #     exit()
    ma_rm = len(dits_rm)/len(dit_qgram_mix)
    ma_add = len(dits_add)/len(dit_qgram_mix)
    ma_mix = len(dits_mix)/len(dit_qgram_mix)
    print(ma_rm,ma_add,ma_mix)

    thor = 0.25
    ave = 0  
    z=0
    while 1:
        count=0
        for i in dits_rm:
            if i>z:
                count+=1
        ma2 = count/len(dit_qgram_mix)
        if ma_rm-ma2>thor:
            print(ma2)
            break
        z += 1
    print(z)

    ave+= z
    z=0
    while 1:
        count=0
        for i in dits_rm:
            if i>z:
                count+=1
        ma2 = count/len(dit_qgram_mix)
        if ma_add-ma2>thor:
            print(ma2)
            break
        z += 1
    print(z)

    ave+= z
    z=0
    while 1:
        count=0
        for i in dits_rm:
            if i>z:
                count+=1
        ma2 = count/len(dit_qgram_mix)
        if ma_mix-ma2>thor:
            print(ma2)
            break
        z += 1
    ave+= z
    print(z)
    # print(begin,z,(begin+z)//2,count_elements_greater_than(dits_rm, 0)/len(dit_qgram_mix),count_elements_greater_than(dits_rm, begin)/len(dit_qgram_mix))
    # exit()
    return ave/3,count_elements_greater_than(dits_rm, (ave/3+z)//2)/len(dit_qgram_mix)

def get_qgram_nums(dit_qgram_change):
    dits = []
    for k in dit_qgram_change:
        if k=='max':
            continue
        l=0
        
        if '#' in k:
            l +=k.count('#')
        if '$' in k:
            l +=k.count('$')
        if l==1 and dit_qgram_change[k]!=0:
            dits.append(dit_qgram_change[k])
    return dits
#产生明文特征
def gen_plain_feature(dit_qgram_mix):
    """
        生成明文的特征,用于后续和密文进行相似度比较,对齐密文和明文
        输入:
            dit_qgram_mix: 明文的变化情况矩阵   [[减少的qgram],[增加的qgram],[减少的qgram中和其他实体减少的qgram相同的变化],[增加的qgram中和其他实体增加的qgram相同的变化],[减少的qgram和其他实体增加的qgram相同的变化]]
        输出:
            dit_palin_feature:      明文的特征,根据dit_qgram_mix[2]\[3]\[4]得来

    """
    dit_palin_feature = {}

    max_f_rm = 0
    for i in dit_qgram_mix:
    
        lst_feature = {}
  
        f_rm_1,f_rm_2,max_qgram_rm = gen_gram_feature(dit_qgram_mix[i][2])
       
        f_add_1,f_add_2,max_qgram_add = gen_gram_feature(dit_qgram_mix[i][3])


        f_mix_1,f_mix_2,max_qgram_mix= gen_gram_feature(dit_qgram_mix[i][4])
   
        dit_palin_feature[i]=[f_rm_1,f_rm_2,f_add_1,f_add_2,f_mix_1,f_mix_2,max_qgram_rm,max_qgram_add,max_qgram_mix]
      
        if i=='BY292032$':
            print(i,dit_palin_feature[i])
        # break
        max_f_rm = max(len(f_rm_1),max_f_rm)
    return dit_palin_feature,max_f_rm



def gen_gram_feature(dit_qgram_change):
    """
    统计以生成明文特征
    输入:
        dit_qgram_change:明文中每个实体和其他实体相同的变化情况
    输出：
        相同变化的频率信息,特征的数量,最多有几个qgram和其他实体同时变化
    """

    f_1 = []
    f_2 = 0

    for k in dit_qgram_change:
        if '#' in k:
            l =k.count('#')
        else:
            l =k.count('$')
        
        if l!=1:
            for i in range(0,len(k),3):
                ai = k[i:i+3]
                if ai in dit_qgram_change:
                    if dit_qgram_change[ai]<global_a:
                        continue
                    dit_qgram_change[ai]+=dit_qgram_change[k]

    max_qgram = 0
    for k in dit_qgram_change:
        if '#' in k:
            l =k.count('#')
        else:
            l =k.count('$')

  

        if max_qgram<l:
            max_qgram =l
        if l>1:
                continue 
  
        if dit_qgram_change[k]>global_a:
            f_2+=1
            f_1.append(dit_qgram_change[k])
  
    f_1.sort(reverse=True)

    return f_1,f_2,max_qgram


#产生密文特征
def gen_enc_feature(table_bit_rm,table_bit_add,table_bit_mix,dit_bit_mix):
    # print("\n",table_bit_rm,'\n')
    lst_feature = {}
    tmp_dit = {}
    for i in dit_bit_mix:

        f_rm_1,f_rm_2,max_bit_rm,a,tmp_bit_dit = gen_feature(dit_bit_mix[i][3],table_bit_rm,dit_bit_mix[i][6])

        f_add_1,f_add_2,max_bit_add,a1,tmp_bit_dit2 = gen_feature(dit_bit_mix[i][4],table_bit_add,dit_bit_mix[i][7])

        f_mix_1,f_mix_2,max_bit_mix,a2,tmp_bit_dit3 = gen_feature(dit_bit_mix[i][5],table_bit_mix,dit_bit_mix[i][8])

        lst_feature[i]=[f_rm_1,f_rm_2,f_add_1,f_add_2,f_mix_1,f_mix_2,max_bit_rm,max_bit_add,max_bit_mix,a,a1,a2]
        tmp_dit[i] = [tmp_bit_dit,tmp_bit_dit2,tmp_bit_dit3]
        
        # # print(i,lst_feature[i],tmp_lst_rm)
        # if i[0:-1] =='AA188158':
        #     print(i,lst_feature[i],tmp_lst_rm)
           
    return lst_feature,tmp_dit

def gen_feature(dit_bit_change,table_bit,tmp_bit,global_a):
    f_1 = []
    f_2 = 0
    max_bit = 0

    tmp_bit_dit = {}
    for k in dit_bit_change:
       
        max_bit=max(dit_bit_change[k][1],max_bit)
        #去除容易被哈希冲突影响的特征
 
        if  dit_bit_change[k][0]<global_a:
            continue

        tmp_bit_dit[k] = [tmp_bit[k][1],dit_bit_change[k][0]]
        if dit_bit_change[k][0]>global_a:
            f_2+=1
        
            f_1.append(dit_bit_change[k][0])
   
    f_1.sort(reverse=True)

    a= max_bit
    if dit_bit_change==0:
        max_bit=[0]
    else:
        if max_bit not in table_bit:
            max_bit=[0]
        else:
            max_bit = table_bit[max_bit]
  

    return f_1,f_2,max_bit,a,tmp_bit_dit




def get_sim_lst(lst1,lst2,n1,n2,max_qgram,max_bit,global_b): 
    if n1!=0 or n2!=0:
        if len(set(lst1).intersection(set(lst2)))/max(n1,n2)>=0.6:
            return 1

    if max_qgram in max_bit:
        sim = 1
    else:
        sim = 1
        
    lsts = []
    
    for i,v in enumerate(lst2):
        max_index = -9
        max_max =9999
        if len(lsts)==len(lst1):
            break
        for j,v2 in enumerate(lst1):
            if j in lsts:
                continue
            if max_max >abs(v2-v):
                max_max = abs(v2-v)
                max_index=j
    
        lsts.append(max_index)
        s = abs(lst1[max_index]-lst2[i])
        sim += max(1-global_b*s/max(lst1[max_index],lst2[i]),0)
    
    d = abs(n1-n2)
    sim += d*0.5+-0.3*(d)*(d-1)/2

    return sim/(max(n1,n2)+1)
 


def get_bit_qgram(dit,dit_plain_same,dit_enc_total):
    bit_qram = {}
    bit_qgram = {}
    for i in dit:
     
        bit_qram = gen_bit_gram(bit_qram,dit[i][0][0],dit[i][1][0])
        bit_qram =gen_bit_gram(bit_qram,dit[i][0][1],dit[i][1][1])
        
    for j in bit_qram:
        tmplst = []
        for k in bit_qram[j][0]:
            if bit_qram[j][0][k]/bit_qram[j][1]>0.8:
                tmplst.append(k)
        bit_qgram[j]=tmplst
    
  
    return bit_qgram

      
    
def gen_bit_gram(bit_qram,qgaram,bit):
    for gram in qgaram:
        c1 = qgaram[gram]
        if c1<10:
            continue
        min_gap = 9999
        min_index=999
        for b in bit:
            c2 = bit[b][1]
            if min_gap>abs(c2-c1):
                min_gap=abs(c2-c1)
                min_index=b
        if gram not in bit_qram:
            t = {}
            for z in bit[min_index][0]:
                t[z]=1

            bit_qram[gram]=[t,1]
            
        else:
            for z in bit[min_index][0]:
                if z not in bit_qram[gram][0]:
                    bit_qram[gram][0][z]=1
                else:
                    bit_qram[gram][0][z]+=1
            bit_qram[gram][1]+=1

    return bit_qram




def get_dit_qgram(result):
    """
    根据已知密文和明文的对应关系识别出不同的qgram可能会被哈希到那些位置
    result 已知的密文对应的明文  [明文1.明文2,密文1,密文2]
    返回一个q-gram和密文bit位位置对应的字典 dit[qgram] = [包含该qgram的明文数据]
    """
    qm1=q-1
    qgram_dit = {}
    
    # 根据已知明文的构建q-gram字典 dit[qgram] = [包含该qgram的明文数据] 2:rm add
    for key in range(1):
        qgram_feature = {}
        tmp_qgram_dit = {}
        for i in result:
            # print(i,result[i])
            # print(dit_bit_mix[result[i][-1]])
            # print(dit_qgram_mix[i+'$'])
            
            # exit()
            atr1 = result[i][key][0]
            atr2 = result[i][key][1]
        
            qgram_a = set([atr1[i:i+q]+attr_symbol[0] for i in range(len(atr1) - qm1)]+[atr2[i:i+q]+attr_symbol[1] for i in range(len(atr2) - qm1)])
            dits = {}
            for j in qgram_a:
                dits[j]=1
            for j in qgram_a:
                if j not in qgram_feature:
                    qgram_feature[j]=[qgram_a,i]
                else:
              
                    qgram_feature[j][0] = qgram_a.intersection(qgram_feature[j][0])
                    qgram_feature[j]+=[i]
    
        for i in qgram_feature:
            tmp = {}
            lst = bitarray.bitarray(l)
            lst.setall(0)
            if len(qgram_feature[i])<=3:
        
                continue
            if len(qgram_feature[i][0])>=2:
                continue
            for j in qgram_feature[i][1:]:
                for k in range(l):

                    if result[j][key+2][k]=='1':
                        if str(k) not in tmp:
                            tmp[str(k)]=1
                        else:
                            tmp[str(k)]+=1

            for j in tmp:
                if tmp[j]/len(qgram_feature[i][1:])>0.9:
                    lst[int(j)]=1
            if lst.count()==0:
                continue
            tmp_qgram_dit[i]=[lst,len(qgram_feature[i][1:]),qgram_feature[i][0]]
        tmp_qgram_dit = dict(sorted(tmp_qgram_dit.items(),key=lambda x:len(x[1][0]),reverse=True))
        tmp_qgram_dit2 = dict(sorted(tmp_qgram_dit.items(),key=lambda x:len(x[1][0]),reverse=True))

        for i in tmp_qgram_dit:
            # print(tmp_qgram_dit[i][0])
            for j in tmp_qgram_dit2:
                if tmp_qgram_dit2[j][0].count()==0:
                    continue
                if i!=j:
                    if (tmp_qgram_dit[i][0] & tmp_qgram_dit2[j][0]).count() /tmp_qgram_dit2[j][0].count()>0.8:
                        tmp_qgram_dit[i][0] = tmp_qgram_dit[i][0] ^ (tmp_qgram_dit[i][0] & tmp_qgram_dit2[j][0])


            # print(tmp_qgram_dit[i][0])
        #     # exit()
        # c= 0 
        # for i in tmp_qgram_dit:
            
        #     if tmp_qgram_dit[i][1]!=0:
        #         # print(list(tmp_qgram_dit[i][0].itersearch(1)))
        #         for j in tmp_qgram_dit[i][2]:

                        
        #             if j not in tmp_qgram_dit:
        #                 continue
        #             if tmp_qgram_dit[j][0].count()==0:
        #                 continue
        #             # if tmp_qgram_dit[j][1]==1:
        #             #     continue
        #             if i!=j:
        #                 if (tmp_qgram_dit[i][0] & tmp_qgram_dit[j][0]).count() /tmp_qgram_dit[j][0].count()>0.8:
        #                     tmp_qgram_dit[i][0] = tmp_qgram_dit[i][0] ^ (tmp_qgram_dit[i][0] & tmp_qgram_dit2[j][0])
                        
                # print(list(tmp_qgram_dit[i][0].itersearch(1)),'\n')
   
        # print(c)
        # exit()
        for i in tmp_qgram_dit:
            if tmp_qgram_dit[i][0].count()==0:
                continue
            if i not in qgram_dit:
                qgram_dit[i]=tmp_qgram_dit[i]
                

    return qgram_dit
def count_palin_qgram(qgram_dit,plain,unmatch_plain):
    dit = {}
    qm1=q-1


    for i in unmatch_plain:
        
        atr1 = unmatch_plain[i][0][0]
        atr2 = unmatch_plain[i][0][1]
        if type(atr1)==float  :
            atr1 = ''

        if type(atr2)==float  :
            atr2 =''
        qgram_a = set([atr1[i:i+q]+attr_symbol[0] for i in range(len(atr1) - qm1)]+[atr2[i:i+q]+attr_symbol[1] for i in range(len(atr2) - qm1)])

        for j in qgram_a:
            dit[j]=1

    print('未识别的密文有多少qgram:',len(dit))   
    for i in qgram_dit:
        dit[i]=1

    print('总共有多少qgram:',len(dit))
    for i in dit_plain_same:
    
        atr1 = dit_plain_same[i][0][0]
        atr2 = dit_plain_same[i][0][1]
        if type(atr1)==float  :
            atr1 = ''

        if type(atr2)==float  :
            atr2 =''
        qgram_a = set([atr1[i:i+q]+attr_symbol[0] for i in range(len(atr1) - qm1)]+[atr2[i:i+q]+attr_symbol[1] for i in range(len(atr2) - qm1)])

        for j in qgram_a:
            dit[j]=1

    dit_qgram = {}
    t = 0
    for j in dit:
        dit_qgram[j]=t
        t+=1

    print('识别出多少qgram',len(qgram_dit))
    return dit_qgram

def identification_enc_data(qgram_dit,dit_plain_same,dit_enc_total,dit1):
    """
        # 根据dit字典算出每个密文可能会包含的qgram，进而和明文进行匹配
        输入：
            dit_enc_part   未知的密文 [密文1 密文2]
            dit_plain_same 没有变化的明文数据  [明文1 , 明文2]
    """
    qm1=q-1
    palin_qgram_feature = {}
    dit_match = {}

    # 将每个明文转化为qgram的形式
    begin = time()
    for i in dit_plain_same:
        
        atr1 = dit_plain_same[i][0][0]
        atr2 = dit_plain_same[i][0][1]
        if type(atr1)==float  :
            atr1 = ''

        if type(atr2)==float  :
            atr2 =''
        qgram_a = set([atr1[i:i+q]+attr_symbol[0] for i in range(len(atr1) - qm1)]+[atr2[i:i+q]+attr_symbol[1] for i in range(len(atr2) - qm1)])

        palin_qgram_feature[i]=qgram_a
 
    # 根据dit字典算出每个密文可能会包含的qgram
    palin = {}
    for i in dit_enc_total:
        if i not in dit_plain_same:
            continue
        sr = bitarray.bitarray(dit_enc_total[i][0])
     
        c2 = []
        for j in qgram_dit:
            if (sr & qgram_dit[j][0]).count()/qgram_dit[j][0].count()>=0.8:
                c2.append(j)
        palin[i] = [c2,dit_plain_same[i][0][0],dit_plain_same[i][0][1],palin_qgram_feature[i]]

    plain_enc = {}

    block = minlsh(palin,palin_qgram_feature)
    begin = time()

    dit_plain_bitarr = {}
    for i in palin_qgram_feature:
        b = bitarray.bitarray(len(dit1))
        b.setall(0)
        for k in palin_qgram_feature[i]:
            b[dit1[k]]=1
        dit_plain_bitarr[i]=b

    one_to_more = {}
    none_match = 0
    for i in palin:

        a = bitarray.bitarray(len(dit1))
        a.setall(0)
        tmp_origin =defaultdict(lambda:[])
        for j in palin[i][0]:
            a[dit1[j]]=1

        if len( block[i])==0:
            none_match+=1
            continue

        for j in block[i]:
            intersection = (a & dit_plain_bitarr[j]).count()
            union = (a | dit_plain_bitarr[j]).count()
            s =  2 * intersection / (union + intersection)
            tmp_origin[j]=s

        tmp_origin = dict(sorted(tmp_origin.items(),key=lambda x:x[1],reverse=True))

        lst = []
        k=0
        tmp = 0
        for j in tmp_origin:
            if k>10:
                break
            lst.append(j)

            if tmp!=tmp_origin[j]:
                k+=1
            tmp=tmp_origin[j]

        one_to_more[i] = lst
        # print(len(lst),m1,palin[i][0],palin_qgram_feature[mindex]) 
        plain_enc[i] = lst[0]
        dit_match[i]=dit_plain_same[lst[0]]+dit_enc_total[i]
        palin[i][0].sort()
        # print(i,palin[i],mindex)
    print('匹配明文的时间',time()-begin,'\n')

    t = 0
    more = 0
    error = 0
    for i in plain_enc:
        if i==plain_enc[i] or palin_qgram_feature[i]==palin_qgram_feature[plain_enc[i]]:
            t+=1
        elif i in one_to_more[i]:
            more+=1
        else:
            error+=1
            # print(i,plain_enc[i],palin[i][0],palin_qgram_feature[i],palin_qgram_feature[plain_enc[i]])

    # print(len(palin),len(plain_enc),len(dit_plain_same))

    return t,dit_match,more,error,none_match
   

def checkgram(qgram_dit,dit_grama):
    c=0
    eror= 0
    z=0
    zz = 0
    zz2=0

    zz3=0
    for i in qgram_dit:
        s=0
        ai = list(qgram_dit[i][0].itersearch(1))
        # print(ai)
        if i in dit_grama:
            # print(dit_grama[i],ai)
            s = len(set(dit_grama[i]).intersection(set(ai)))/len(set(ai))
            zz+=s
            zz2+=1
            zz3+=len(set(dit_grama[i]).intersection(set(ai)))/len(set(dit_grama[i]))
            # if s/len(set(ai))<0.6:
            #     s=0
            # print(i,set(dit_grama[i]),set(ai))
        # print(s)
        if s>0.6:
            c+=1
        else:
            eror+=1
        if i not in dit_grama:
            z+=1
        # else:
        #     if i in dit_grama:
        #         print(i,set(dit_grama[i]),set(qgram_dit[i][0]))
        #     elif i in dit_gramb:
        #         print(i,set(dit_gramb[i]),set(qgram_dit[i][0]))
    print(zz/zz2,zz3/zz2)
    print('准确率',c,len(qgram_dit),c/len(qgram_dit),eror,z)

def get_same_qgram(result,dit_plain_same,dit_enc_total,unmatch_plain,unmatch_bit,total_change,dit_grama,p):
    """
    result 已知的密文对应的明文  [明文1.明文2,密文1,密文2]
    dit_enc_part   未知的密文 [密文1 密文2]
    dit_plain_same 没有变化的明文数据  [明文1 , 明文2]
    """
    # 得到qgram对应的bit位置
    qgram_dit = get_dit_qgram(result)
   
    match_nums=0


    dit1 = count_palin_qgram(qgram_dit,dit_plain_same,unmatch_plain)

    checkgram(qgram_dit,dit_grama)

    match_nums2,dit_match,cor_1_to_m_un,error2,none_match = identification_enc_data(qgram_dit,unmatch_plain,unmatch_bit,dit1)

    qgram_dit2 = get_dit_qgram(dit_match)
    for i in qgram_dit2:
        if i not in qgram_dit:
            qgram_dit[i]=qgram_dit2[i]

    f = open('this-qgram_dit-euro-tmh.txt', 'a',encoding='gbk')      
    for i in qgram_dit:
        s=''
        indexes = str(qgram_dit[i][0].to01())
        # 将字符串写入到文件中
        f.write(i+'-'+indexes+'\n')
        
    # 关闭文件
    f.close()

    dit1 = count_palin_qgram(qgram_dit,dit_plain_same,unmatch_plain)
    checkgram(qgram_dit,dit_grama)
    
    if len(dit_plain_same)!=0:
        
        match_nums,dit_match2,cor_1_to_m,error,none_match2 = identification_enc_data(qgram_dit,dit_plain_same,dit_enc_total,dit1)
    # print(match_nums,total_nums,pre,nums1,pre2)

    s =len(result)

    
    print(f'变化记录总数{len(unmatch_plain)+s},变化记录中识别的总数：{s},其中正确的：{p}，错误{s-p}  \
          其中未识别的的：{len(unmatch_plain)},识别后正确的 1-1:{match_nums2},  1-m:{cor_1_to_m_un} 错误{error2}  未进行匹配的{none_match}')

    print(f'没有变化的明文数量{len(dit_plain_same)},正确匹配的{match_nums},一对多匹配{cor_1_to_m},错误的数目{error} 未进行匹配的{none_match2}')

    print(f' 综合匹配过程:\
            总的一对多的数目{cor_1_to_m+cor_1_to_m_un}\
            总的1-1匹配数目:{match_nums2+match_nums+p},\
            错误匹配数目：{error2+error+s-p}\
            未匹配的数目：{none_match+none_match2}\
            总的数目是:{len(dit_enc_total)+len(unmatch_bit)+s},\
            1-1准确率：{(match_nums2+match_nums+p)/(len(dit_plain_same)+len(unmatch_plain)+s)}\
            一对多的准确率:{(cor_1_to_m+cor_1_to_m_un)/(len(dit_plain_same)+len(unmatch_plain)+s)}')
    
def get_unmatch_record(result,dit_plain_change,dit_enc_part):
    """
        由于舍弃了一部分召回率在此补正
        得到没有被匹配的部分变化数据
        result:已知匹配的密文对应的明文  [明文1.明文2,密文1,密文2]
        dit_plain_change: 所有的明文
        dit_enc_part: 所有的密文

        返回
        未匹配的明文，未匹配的密文
    """
    unmatch_plain ={}
    unmatch_bit ={}
    for i in dit_plain_change:
        if i not in result:
            unmatch_plain[i]=dit_plain_change[i]
            if i in dit_enc_part:
                unmatch_bit[i] = dit_enc_part[i]
   
    return unmatch_plain,unmatch_bit

def stable_marriage(men_pref, women_pref):
    """
    使用Gale-Shapley算法解决稳定婚姻匹配问题

    参数:
        men_pref (dict): 男方节点的偏好字典，键为男方节点，值为女方节点列表，按偏好从高到低排序
        women_pref (dict): 女方节点的偏好字典，键为女方节点，值为男方节点列表，按偏好从高到低排序

    返回:
        dict: 匹配结果，键为男方节点，值为与其匹配的女方节点
    """
    matches = {}  # 存储匹配结果的字典
    unmatched_men = list(men_pref.keys())  # 未匹配的男方节点列表

    while unmatched_men:
        man = unmatched_men[0]  # 取出一个未匹配的男方节点

        # 获取男方节点的偏好列表和当前偏好索引
        pref_list = men_pref[man]
        pref_index = len(matches[man]) if man in matches else 0

        # 获取当前男方节点在偏好列表中的下一个偏好女方节点
        woman = pref_list[pref_index]

        # 如果女方节点未匹配，则将男方节点与女方节点匹配
        if woman not in matches:
            matches[man] = woman
            unmatched_men.pop(0)  # 将已匹配的男方节点从未匹配列表中移除
        else:
            # 如果女方节点已匹配，比较当前男方节点和原匹配男方节点在女方节点的偏好中的位置
            curr_man = matches[woman]
            curr_pref_index = women_pref[woman].index(curr_man)
            new_pref_index = women_pref[woman].index(man)

            # 如果当前男方节点在女方节点的偏好中的位置更高，则保持原匹配，否则将男方节点与女方节点匹配
            if new_pref_index < curr_pref_index:
                matches[man] = woman
                unmatched_men.pop(0)  # 将已匹配的男方节点从未匹配列表中移除

    return matches

def batcher_iter(data, batch_size):
    """
    批处理迭代器，将数据划分为指定的 batch_size 批次
    """
    for i in range(0, len(data), batch_size):
        yield data[i:i+batch_size]

def compare_bit_gram_feature(dit_palin_feature,dit_enc_feature):
    begin_time     =   time()
    blocks =compare_minlsh(dit_enc_feature,dit_palin_feature)
    print('lsh哈希时间',time()-begin_time)

    keys_list = list(dit_palin_feature.keys())
    begin_time     =   time()   
    epoch = epochs
    batch_size = int(len(dit_palin_feature)/epoch)
    # compare_fearture(keys_list,dit_enc_feature,dit_palin_feature,global_a,global_b,blocks)
    with multiprocessing.Pool(epoch) as  pool:
        for batch in batcher_iter(keys_list, batch_size):
            this_result =  pool.apply_async(compare_fearture,(batch,dit_enc_feature,dit_palin_feature,global_a,global_b,blocks),callback=save_compre_result)
    
        pool.close()
        pool.join()

    end_time =   time()
    print('二部图匹配 计算相似度的时间',end_time-begin_time)


    tt = {}
    for i in block_smh_result:
        if block_smh_result[i][0] not in tt:
            tt[block_smh_result[i][0]]=[i,block_smh_result[i][1]]
        else:
            if block_smh_result[i][1]>tt[block_smh_result[i][0]][1]:
                tt[block_smh_result[i][0]]=[i,block_smh_result[i][1]]
    result = {}
    for i in tt:
        result[tt[i][0][0:-1]]=dit_plain_change[tt[i][0][0:-1]]+dit_enc_change[i[0:-1]]+[i]

    k=0
    for i in tt:
        if i[:-1]==tt[i][0][0:-1]:
            k+=1
    # block =block_result
    # dit_a = {}
    # for i in block:
    #     for j in block[i]:
    #         if j not in dit_a:
    #             dit_a[j] = {i:block[i][j]}
    #         else:
    #             dit_a[j][i]=block[i][j]
    
    # x1 = {}
    # x2 = {}
    # for i in block:
    #     block[i] = dict(sorted(block[i].items(),key=lambda x:x[1],reverse=True))
    #     x1[i] = list(block[i].keys())
       
    # for i in dit_a:
    #     dit_a[i] = dict(sorted(dit_a[i].items(),key=lambda x:x[1],reverse=True))
    #     x2[i] = list(dit_a[i].keys()) 
      
    # k=0
    # for i in x1:
    #     if i[0:-1]+'#' in x2:
    #         if x1[i][0][0:-1]==i[0:-1] and  x2[i[0:-1]+'#'][0][0:-1]==i[0:-1]:
    #             k+=1
    # print(len(x1),k)

   
    # # 调用stable_marriage函数解决稳定婚姻匹配问题
    # matching = stable_match(x1, x2)
   
    # result = {}
    # gg = 0
    # gg2=0
    # # 输出匹配结果
    # for man, woman in matching:
     
    #     result[man[0:-1]]=dit_plain_change[man[0:-1]]+dit_enc_change[woman[0:-1]]
        
    #     gg2+=1
    #     if man[0:-1]==woman[0:-1]:
    #         gg+=1
    #     # else:
    #     #     print(man,woman,dit_palin_feature[man],dit_enc_feature[woman],dit_enc_feature[man[0:-1]+'#'])
    #     # print(f"{man} is matched with {woman}")
    # print(gg2,gg)
    print(len(tt),k,len(dit_palin_feature))
    
    unmatch_plain,unmatch_bit = get_unmatch_record(result,dit_plain_change,dit_enc_change)
    # #已知变化的密文对应的明文 ，提取有效信息继而识别全部密文
   
    return result,unmatch_plain,unmatch_bit,k

def save_compre_result(result):

    for i in result[0]:
        block_result[i] = result[0][i]
    for i in result[1]:
        block_smh_result[i] = result[1][i]
    return 1


def check_sim(lst1,lst2,global_a):
    if lst1[1]>0 and lst2[1]>0 and lst1[3]>0 and lst2[3]>0 and lst2[5]>0 and lst1[5]>0:
        k=0
        a=0

        if abs(lst1[0][a]-lst2[0][a])<global_a/2:
            k+=1
        if abs(lst1[2][a]-lst2[2][a])<global_a/2:
            k+=1
        if  abs(lst1[4][a]-lst2[4][a])<global_a/2:
            k+=1
        if k>=1:
            return 1
    return 0

def compare_fearture(batch,dit_enc_feature,dit_palin_feature,global_a,global_b,blocks):
    block = {}
    block_smh = {}
    a_th=0.75
    for i in batch:
        index= {} 

        ks=0
        if i not in blocks:
            continue
        for k in blocks[i]:
  
            # if abs(dit_palin_feature[i][1]-dit_enc_feature[k][1])+abs(dit_palin_feature[i][3]-dit_enc_feature[k][3])+abs(dit_palin_feature[i][5]-dit_enc_feature[k][5])>3:
            #     continue

            ks+=1
            s1_1 = max(1,dit_palin_feature[i][1])
            s2_1 = max(1,dit_palin_feature[i][3])
            s3_1 = max(1,dit_palin_feature[i][5])

            s1 = get_sim_lst(dit_palin_feature[i][0],dit_enc_feature[k][0],dit_palin_feature[i][1],dit_enc_feature[k][1],dit_palin_feature[i][6],dit_enc_feature[k][6],global_b)
            if s1<((a_th)*(s3_1+s2_1+s1_1)-s2_1-s3_1)/s1_1:
                continue
           
            s2 = get_sim_lst(dit_palin_feature[i][2],dit_enc_feature[k][2],dit_palin_feature[i][3],dit_enc_feature[k][3],dit_palin_feature[i][7],dit_enc_feature[k][7],global_b)
            if s2<((a_th)*(s3_1+s2_1+s1_1)-s1_1*s1-s3_1)/s2_1:
                continue      
                    
            s3 = get_sim_lst(dit_palin_feature[i][4],dit_enc_feature[k][4],dit_palin_feature[i][5],dit_enc_feature[k][5],dit_palin_feature[i][8],dit_enc_feature[k][8],global_b)
            if s3<0.3:
                continue             
            s  = (s1*s1_1+s2*s2_1+s3*s3_1)/(s3_1+s2_1+s1_1)
            if s>a_th:
                index[k] = s
    
        index = dict(sorted(index.items(),key=lambda x:x[1],reverse=True))


        maxs= 0
        for k in index:
            maxs=k
            if index[k]>a_th:
                block_smh[i]=[k,index[k]]
            break
        if maxs not in index:
            continue
        if index[maxs]>a_th:
            block[i]=index


    return [block,block_smh]
   
def stable_match(men, women):
    free_men = deque(men)
    engaged = defaultdict(lambda: None)
    while free_men:
        i = free_men.popleft()
        # man proposes women according his preferences
        for j in men[i]:
            preference = women[j].index
            fiance = engaged[j]
            # woman accepts the better offer
            if not fiance or preference(i) < preference(fiance):
                engaged[j] = i
                fiance and free_men.append(fiance)
                break
    return [(m, w) for w, m in engaged.items()]

def get_gram_dit(qgram_dit):
    dit_grama ={}
    for i in qgram_dit:
        i[1] = (i[1][1:-1]).split(',')
        for j,v in enumerate(i[1]):
            i[1][j] = int(v)

        dit_grama[i[0]] = i[1]
    
    return dit_grama


def Gen_seq_feature(dit_enc_feature):
    lst_feature2 = {}
    for i in dit_enc_feature[0]:
        if i not in dit_enc_feature[0]:
            continue
        if i not in dit_enc_feature[1]:
            continue
        if i not in dit_enc_feature[2]:
            continue        
        rm = dit_enc_feature[0][i]
        add = dit_enc_feature[1][i]
        mix = dit_enc_feature[2][i]
        lst_feature2[i] = [rm[0],rm[1],add[0],add[1],mix[0],mix[1],rm[2],add[2],mix[2]]
    
    return lst_feature2
    
if __name__ ==   '__main__':

    #==================================base：参数设置===============================================

    path1 = "2/euro-a-rbf-.csv"
    path2 = "2/euro-b-rbf-.csv"
    path_qgram = "2/qgram_dit-euro-rbf.csv"             # qgram字典的路径，用于计算精度等
    attr_index  = [1,2]                                 # 使用的属性 两个 firstname lastname
    q = 2                                               # q-gram
    l = 1000                                            # Length of Bloom Filter
    epochs = 12                                         # Utilizing multithreading to speed up the runtime, for debugging purposes only. 
    thor_signature = 0.7                                # Minimum threshold for signatures added to clusters when adjusting clustering

    if 'rbf' in path1 or 'salt' in path1:               # The location of the hash of the q-grams is different for the different attributes of rbf and salting.
        attr_symbol = ['$','#']
    else:
        attr_symbol = ['$','$']

    #==================================每个qgram对应的bf位位置，用于计算精度===========================
        
    qgram_dit,c_name    = loadCSVData(path_qgram)   # 读取字典内容转换列表
    rec_list,c_name     = loadCSVData(path1)        # 读取csv文件内容转换列表
    rec_list2,c_name    = loadCSVData(path2)        # 读取csv文件内容转换列表

    dit_grama = get_gram_dit(qgram_dit)             # 字典列表 dit['qgram']=[position of bits]
    
    #==================================step1:对记录进行分类===============================================
    times=time()

    dit_plain_same,dit_plain_change,dit_plain_none       = gen_dict_plain([rec_list,rec_list2])      # 没有改变的明文记录  发生改变的记录 不是同一个实体的记录
    dit_enc_same,dit_enc_change,dit_enc_none             = gen_dict_encode([rec_list,rec_list2])     # 没有改变的密文记录  发生改变的密文 不是同一个实体的记录
    
    print('step1 time:',time()-times,'\n')

    #==================================step2:提取变化数据的特征==============================================

    times=time()

    dit_bit_add,dit_bit_rm,dit_bit_mix                 = get_enc_bit(dit_enc_change)                    # 密文, 只增加bit位，只减少bit位和混合
    dit_qgram_add,dit_qgram_rm,dit_qgram_mix           = get_plain_qgram(dit_plain_change)              # 明文，只增加qgram，只减少q-gram和混合

    print('step2 time:',time()-times,'\n')

    
    #==================================step3-1:产生qgram的图==============================================
    times=time()
    
    len_qgram_change    = {'rm':[],'add':[],'mix':[]}
    conter_lst          = {'rm':[],'add':[],'mix':[]}

    gen_qgram_grap(dit_qgram_mix)

    print('step3-1 time:',time()-times,'\n')
   
    #==================================step3-2:产生一个全局的key，舍去掉一些不明显的特征==============================================
    times=time()
 
    global_a,global_b = get_gobal_key(dit_qgram_mix)
    # global_a =120                                     # Can be set based on experience
    print('使用的全局key:',global_a,'. 单个qgram出现的频率小于此值将被丢弃','\n')
    print('step3-2 time:',time()-times,'\n')

    #==================================step3-3:产生qgram的特征==============================================
    times=time()
    max_f_rm = 0                 # Limit the number of clusters to speed up, optional.
    dit_palin_feature,max_f_rm  = gen_plain_feature(dit_qgram_mix)
    print('step3-3 max_f_rm:',max_f_rm)
    print('step3-3 time:',time()-times,'\n')
   

    #==================================step4:产生bit的图==============================================

    dit_enc_feature = [{} for i in range(3)]
    drop_key_len    = {'rm':[],'add':[],'mix':[]}
    rel_id          = {'rm':{},'add':{},'mix':{}}
    begin_time      = time()

    gen_bit_grap(dit_bit_mix,len_qgram_change['rm'],len_qgram_change['add'],len_qgram_change['mix'])

    end_time =   time()
    print('构建密文特征图的时间:',end_time-begin_time,'\n')
    lst_feature2 = {}
    
    #==================================step4-2:产生bit的特征==============================================
    dit_enc_feature = Gen_seq_feature(dit_enc_feature)
   
    #==================================step5:二部图匹配==============================================

    block_result = {}
    block_smh_result = {}
    result,unmatch_plain,unmatch_bit,k = compare_bit_gram_feature(dit_palin_feature,dit_enc_feature)


    
    #==================================step6:识别q-gram以及对未识别的编码记录进行识别==============================================

    get_same_qgram(result,dit_plain_same,dit_enc_same,unmatch_plain,unmatch_bit,len(dit_plain_change),dit_grama,k)