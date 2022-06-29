import Base.push!, Base.haskey

line_integers(s) = [ parse(Int,x) for x in split(s,r" ") if x!="" ]

forbidden_words() = [line_integers(s)[2:end] for s in readlines()]
forbidden_words(filename) = [line_integers(s)[2:end] for s in readlines(filename)]

add_reverses(l :: Array{Array,1}) = vcat(l,map(reverse,l))

mutable struct Trie
        children :: Array{Array{Int,1},1}
        final :: Array{Bool,1}
        data :: Array{Any,1}
        size :: Int
        arity :: Int
end

function allocate!(tr :: Trie, data = Nothing)
        for c in tr.children; push!(c,0); end
        push!(tr.data, data)
        push!(tr.final, false)
        tr.size += 1
end

function makeTrie(; data = Nothing, arity=2)
        tr = Trie([
                Array{Int}(UndefInitializer(),0)
                for k in 1:arity
        ],
        Array{Bool}(UndefInitializer(),0)
        ,[], 0, arity)
        allocate!(tr,data)
        tr
end

function walk(tr :: Trie, w :: Array{Int,1}; create :: Bool = false, data = Nothing, copydata = false, finaldata = data, finalise = false, start = 1)
        p = start
        created = 0
        for c in w
                np = tr.children[c][p]
                if np == 0
                        if create
                                np = (tr.children[c][p] = allocate!(tr,if copydata; copy(data) else data end))
                                created = np
                        else
                                return 0
                        end
                end
                p = np
        end
        if created == p
                tr.data[p] = if copydata; copy(finaldata) else finaldata; end
        end
        if finalise; tr.final[p] = true; end
        p
end

function push!(tr :: Trie, w :: Array{Int, 1}; data = Nothing, finaldata = data, copydata = false)
        n = walk(tr, w, create=true, data=data, copydata = copydata, finaldata = finaldata)
        tr.final[n] = true
        tr
end

function haskey(tr :: Trie, w :: Array{Int,1})
        return walk(tr, w) != 0
end

function digital_words(n; arity=2)
        if n == 0
                return [Array{Int,1}(UndefInitializer(), 0)]
        end
        res = Array{Array{Int,1},1}(UndefInitializer(), 0)
        for w in digital_words(n-1,arity=arity)
                for c in 1:arity
                        push!(res, vcat(w,[c]))
                end
        end
        res
end

function permitted_words(words, forbidden; arity=2)
        tr = makeTrie(arity=arity)
        for w in forbidden
                push!(tr, w)
        end
        res = Array{Array{Int,1},1}(UndefInitializer(),0)
        for w in words
                ps = Set{Int}()
                ok = true
                for c in w
                        push!(ps,1)
                        nps = Set{Int}()
                        for p in ps
                                np = tr.children[c][p]
                                push!(nps,np)
                        end
                        delete!(nps,0)
                        for np in nps
                                if tr.final[np]
                                        ok = false
                                end
                        end
                        ps = nps
                end
                if ok
                        push!(res,w)
                end
        end
        res
end

function trie_words(tr :: Trie; data = false, pos = false, start = 1)
        stack = [(start,Int[])]
        res = []
        while length(stack)>0
                p, w = pop!(stack)
                if tr.final[p]
                        if data || pos
                                entry = Any[w]
                                if data
                                        push!(entry, tr.data[p])
                                end
                                if pos
                                        push!(entry, p)
                                end
                                push!(res,(entry...,))
                        else
                                push!(res,w)
                        end
                end
                for c in 1:tr.arity
                        np = tr.children[c][p] 
                        if np != 0
                                push!(stack, (np, vcat(w,Int[c])))
                        end
                end
        end
        res
end

function prefix_tree(words :: Array{Array{Int,1},1})
        arity = max(vcat(words...)...)
        tr = makeTrie(data=[],arity=arity)
        for kk in 1:length(words)
                w = words[kk]
                for k in 1:length(w)
                        n = walk(tr,w[k:-1:1],create=true,data=[], copydata = true, finalise = true)
                        push!(tr.data[n], (kk,k))
                end
        end
        tr
end

suffix_tree(words) = prefix_tree(map(reverse, words))

function encode_affixes(str :: Trie, ptr :: Trie, words :: Array{Array{Int,1},1})
        ls = map(length, words)
        idx = Dict{Tuple,Int}()
        n = 0
        pw = trie_words(ptr, pos=true, data=true)
        sw = trie_words(str, pos=true, data=true)
        for (w,ps,loc) in pw
                new = Array{Int,1}()
                for (wn,l) in ps
                        if l < ls[wn]
                                idx[(wn,ls[wn]-l)] = (n+=1)
                                push!(new, n)
                        end
                end
                ptr.data[loc] = new
        end
        for (w,ps,loc) in sw
                new = Array{Int,1}()
                for (wn,l) in ps
                        if l < ls[wn]
                                push!(new,idx[(wn,l)])
                        end
                end
                str.data[loc] = new
        end
        n
end

function preprocess_forbidden_words(words :: Array{Array{Int,1},1})
        ll = max(map(length,words)...)
        arity = max(vcat(words...)...)
        all_words = digital_words(ll-1, arity = arity)
        good_words = permitted_words(all_words,words, arity = arity)
        str = suffix_tree(words)
        ptr = prefix_tree(words)
        n = encode_affixes(str, ptr, words)
        (words, good_words, ptr, str, n)
end

preprocess_forbidden_words() = preprocess_forbidden_words(forbidden_words())
preprocess_forbidden_words(fname :: String) = preprocess_forbidden_words(forbidden_words(fname))

function word_risks(w :: Array{Int,1}, ptr :: Trie, str :: Trie, n)
        end_prefixes = ones( Int, n)
        start_suffixes = ones(Int, n)
        p = 1
        for c in w
                np = str.children[c][p]
                if np == 0
                        break
                end
                for i in str.data[np]
                        start_suffixes[i] = 2
                end
                p = np
        end
        p = 1 
        for c in reverse(w)
                np = ptr.children[c][p]
                if np == 0
                        break
                end
                for i in ptr.data[np]
                        end_prefixes[i] = 2
                end
                p = np
        end
        (start_suffixes, end_prefixes)
end

function all_risks(ws :: Array{Array{Int,1},1}, ptr :: Trie, str :: Trie, n)
        prtr = makeTrie(data=Int[])
        srtr = makeTrie(data=Int[])
        for wi in 1:length(ws)
                w = ws[wi]
                (ss, ep) = word_risks(w, ptr, str, n)
                pp = walk(prtr, ep, create=true, data=Int[], copydata = true, finalise = true)
                push!(prtr.data[pp], wi)
                sp = walk(srtr, ss, create=true, data=Int[], copydata = true, finalise = true)
                push!(srtr.data[sp], wi)
        end
        (prtr, srtr)
end

struct TrieScaffolding
        subtree_sizes :: Array{Int,1}
        parent :: Array{Int,1}
        offset :: Array{Int,1}
        last_letter :: Array{Int,1}
        depth :: Array{Int, 1}
end

function trie_scaffolding(tr :: Trie)
        n = tr.size
        res = TrieScaffolding(
                zeros(Int,n),
                zeros(Int,n),
                zeros(Int,n),
                zeros(Int,n),
                zeros(Int,n),
        )
        for k in 1:n
                for c in 1:tr.arity
                        p = tr.children[c][k]
                        if p != 0
                                res.parent[p] = k
                                res.last_letter[p] = c
                                res.depth[p] = res.depth[k]+1
                        end
                end
        end
        for k in n:-1:1
                if tr.final[k]
                        res.subtree_sizes[k] += 1
                end
                p = res.parent[k]
                if p != 0
                        res.subtree_sizes[p] += res.subtree_sizes[k]
                end
        end
        for k in 1:n
                s = res.offset[k]
                for c in 1:tr.arity
                        p = tr.children[c][k]
                        if p != 0
                                res.offset[p] = s
                                s += res.subtree_sizes[p]
                        end
                end
        end
        res
end

function trie_intersections(prtr :: Trie, srtr :: Trie)
         psc = trie_scaffolding(prtr)
         ssc = trie_scaffolding(srtr)

         res = makeTrie(data=Int[])
         stack = [(1,Int[],0,1,0,1,1)]

         while length(stack)>0
                (pp,w,ppd,sp,spd,sm,p) = pop!(stack)
                final = true
                if srtr.final[sp]
                        np = walk(res, Int[1], create = true, data = Int[], start = p, copydata = true)
                        nsp = ssc.parent[sp]
                        if nsp != 0
                                final = false
                                push!(stack, (pp, w, ppd, nsp, spd-1, ssc.last_letter[sp]+1,np))
                        end
                elseif spd < ppd
                        found = false
                        for c in sm:2
                                nsp = srtr.children[c][sp]
                                if nsp != 0
                                        if c == w[spd+1] == 2
                                                nrw = ones(Int,ssc.subtree_sizes[nsp]) .* 2
                                                np = walk(res,nrw,create=true,data=Int[],start=p, copydata = true)
                                                push!(stack, (pp, w, ppd, sp, spd, c+1, np))
                                        else
                                                push!(stack, (pp, w, ppd, nsp, spd+1, 1, p))
                                        end
                                        final = false
                                        found = true
                                        break
                                end
                        end
                        if (!found)
                                nsp = ssc.parent[sp]
                                if nsp != 0
                                        final = false
                                        push!(stack, (pp, w, ppd, nsp, spd-1, ssc.last_letter[sp]+1,p))
                                end
                        end
                else
                        found = false
                        final = false
                        for c in 1:2
                                npp = prtr.children[c][pp]
                                if npp!=0
                                        push!(stack, (npp, vcat(w,Int[c]),ppd+1,sp,spd,sm,p))
                                        found = true
                                end
                        end
                end
                if final
                        res.final[p] = true
                        push!(res.data[p],pp)
                end
         end

         res
end

function trie_word(p :: Int, tr :: Trie, scaf :: TrieScaffolding)
        res = Int[]
        while scaf.parent[p] != 0
                push!(res, scaf.last_letter[p])
                p = scaf.parent[p]
        end
        reverse(res)
end

function trie_first_leaf(p :: Int, tr :: Trie)
        while ! tr.final[p]
                for c in tr.children
                        if c[p] != 0
                                p = c[p]
                                break
                        end
                end
        end
        p
end

function compute_compatibility(words :: Array{Array{Int,1},1}; uniquewords = Nothing, multiplicities = Nothing, comultiplicities = Nothing)
        (good_words, ptr, str, splits) = preprocess_forbidden_words(words)[2:end]
        ptrs = trie_scaffolding(ptr)
        strs = trie_scaffolding(str)
        (prtr, srtr) = all_risks(good_words, ptr, str, splits)
        pinter = trie_intersections(prtr,srtr)
        sinter = trie_intersections(srtr,prtr)
        pfcodes = [data[end] for (w,data) in trie_words(pinter,data=true)]
        sfcodes = [data[end] for (w,data) in trie_words(sinter,data=true)]
        pfwords = [prtr.data[trie_first_leaf(k,prtr)][1] for k in pfcodes]
        sfwords = [srtr.data[trie_first_leaf(k,srtr)][1] for k in sfcodes]
        pidx = Dict{Int,Array{Int,1}}()
        sidx = Dict{Int,Array{Int,1}}()
        for (w,data) in trie_words(pinter,data=true)
                all=[]
                for code in data
                        leaves = trie_words(prtr,data=true,start=code)
                        for l in leaves
                                append!(all,l[2])
                        end
                end
                sort!(all)
                pidx[all[1]]=all
        end
        for (w,data) in trie_words(sinter,data=true)
                all = []
                for code in data
                        leaves = trie_words(srtr,data=true,start=code)
                        for l in leaves
                                append!(all,l[2])
                        end
                end
                sort!(all)
                sidx[all[1]]=all
        end
        pcnt=Dict([a=>length(b) for (a,b) in pidx])
        scnt=Dict([a=>length(b) for (a,b) in sidx])
        puw = sort([(good_words[wi],n) for (wi,n) in pcnt])
        suw = sort([(good_words[wi],n) for (wi,n) in scnt])
        if uniquewords != Nothing
                write(uniquewords,
                join([join(w, " ")*"        "*string(n) for (w,n) in suw],"\n")*"\n")
        end
        if multiplicities != Nothing
                write(multiplicities,
                join([string(length(l),pad=8)*"    "*join(l," ") for (wi, l) in sort(collect(sidx))],"\n")*"\n")
        end
        if comultiplicities != Nothing
                write(comultiplicities,
                join([string(length(l),pad=8)*"    "*join(l," ") for (wi, l) in sort(collect(pidx))],"\n")*"\n")
        end
        sort(collect(suw))
end
