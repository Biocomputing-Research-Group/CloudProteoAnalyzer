user_data_path = "the folder path for storing the user configuration files and user data"
host_inf = {"host":"supercomputer address",
                "username":"supercomputer user name",
                "passwords":"supercomputer user password"}
smtp_server = "smtp.gmail.com"
port = 587  
sender_email = "email address for send back results"
password = "email password"

ptm_list = ["PTM{~} = M","PTM{~} = MPK","PTM{!} = NQR","PTM{@} = STYHD","PTM{>to1} = STYHD","PTM{<to2} = ST","PTM{%} = K","PTM{^} = KRED","PTM{&} = KR","PTM{*} = K","PTM{(} = C","PTM{)} = Y","PTM{/} = C","PTM{$} = D"]

ptm_expanation ={"PTM symbol":["PTM{~} = M","PTM{!} = NQR","PTM{@} = STYHD","PTM{>to1} = STYHD","PTM{<to2} = ST","PTM{%} = K","PTM{^} = KRED","PTM{&} = KR","PTM{*} = K","PTM{(} = C","PTM{)} = Y","PTM{/} = C","PTM{$} = D"],"Detail of PTM":["Oxidation of Met","Deamidation of NQ","Phosphorylation","Phosphorylation with losing HPO3","Phosphorylation with losing HPO3 and H2O","Acetylation","Mono-methylation","Di-methylation","Tri-methylation","S-Nitrosylation, search with natural Cys","Nitration","IAA blocking","beta-methylthiolation"]}

enzyme_dict={"Trypsin":["KR","ACDEFGHIJKLMNQRSTVWY"],
                "Trypsin/P":["KR","ACDEFGHIJKLMNPQRSTVWY"],
                "Lys_C":["K","ACDEFGHIJKLMNQRSTVWY"],
                "Lys_N":["ACDEFGHIJKLMNPQRSTVWY","K"],
                "Arg_C":["R","ACDEFGHIJKLMNQRSTVWY"],
                "Asp_N":["ACDEFGHIJKLMNPQRSTVWY","D"],
                "CNBr":["M","ACDEFGHIJKLMNPQRSTVWY"],
                "Glu_C":["DE","ACDEFGHIJKLMNQRSTVWY"],
                "PepsinA":["PL","ACDEFGHIJKLMNQRSTVWY"],
                "Chymotrypsin":["FWYL","ACDEFGHIJKLMNQRSTVWY"]
                }

N_iso_prorata="<ISOTOPOLOGUE name = \"reference\" >\n				<!--	Name	C	H	O	N	P	S	C*	H*	O*	N*	P*	S*	-->\n				<R>	NTerm,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n				<R>	CTerm,	0,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n				<R>	L,	6,	11,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	A,	3,	5,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	S,	3,	5,	2,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	G,	2,	3,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	V,	5,	9,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	E,	5,	7,	3,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	K,	6,	12,	1,	0,	0,	0,	0,	0,	0,	2,	0,	0	</R>\n				<R>	I,	6,	11,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	T,	4,	7,	2,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	D,	4,	5,	3,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	R,	6,	12,	1,	0,	0,	0,	0,	0,	0,	4,	0,	0	</R>\n				<R>	P,	5,	7,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	N,	4,	6,	2,	0,	0,	0,	0,	0,	0,	2,	0,	0	</R>\n				<R>	F,	9,	9,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	Q,	5,	8,	2,	0,	0,	0,	0,	0,	0,	2,	0,	0	</R>\n				<R>	Y,	9,	9,	2,	0,	0,	0,	0,	0,	0,	1,	0,	0	</R>\n				<R>	M,	5,	9,	1,	0,	0,	1,	0,	0,	0,	1,	0,	0	</R>\n				<R>	H,	6,	7,	1,	0,	0,	0,	0,	0,	0,	3,	0,	0	</R>\n				<R>	C,	3,	5,	1,	0,	0,	1,	0,	0,	0,	1,	0,	0	</R>\n				<R>	W,	11,	10,	1,	0,	0,	0,	0,	0,	0,	2,	0,	0	</R>\n"

C_iso_prorata="			<ISOTOPOLOGUE name = \"reference\" >\n		<!--	Name	C	H	O	N	P	S	C*	H*	O*	N*	P*	S*	-->\n				<R>	NTerm,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n				<R>	CTerm,	0,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0	</R>\n				<R>	L,	0,	11,	1,	1,	0,	0,	6,	0,	0,	0,	0,	0	</R>\n				<R>	A,	0,	5,	1,	1,	0,	0,	3,	0,	0,	0,	0,	0	</R>\n				<R>	S,	0,	5,	2,	1,	0,	0,	3,	0,	0,	0,	0,	0	</R>\n				<R>	G,	0,	3,	1,	1,	0,	0,	2,	0,	0,	0,	0,	0	</R>\n				<R>	V,	0,	9,	1,	1,	0,	0,	5,	0,	0,	0,	0,	0	</R>\n				<R>	E,	0,	7,	3,	1,	0,	0,	5,	0,	0,	0,	0,	0	</R>\n				<R>	K,	0,	12,	1,	2,	0,	0,	6,	0,	0,	0,	0,	0	</R>\n				<R>	I,	0,	11,	1,	1,	0,	0,	6,	0,	0,	0,	0,	0	</R>\n				<R>	T,	0,	7,	2,	1,	0,	0,	4,	0,	0,	0,	0,	0	</R>\n				<R>	D,	0,	5,	3,	1,	0,	0,	4,	0,	0,	0,	0,	0	</R>\n				<R>	R,	0,	12,	1,	4,	0,	0,	6,	0,	0,	0,	0,	0	</R>\n				<R>	P,	0,	7,	1,	1,	0,	0,	5,	0,	0,	0,	0,	0	</R>\n				<R>	N,	0,	6,	2,	2,	0,	0,	4,	0,	0,	0,	0,	0	</R>\n				<R>	F,	0,	9,	1,	1,	0,	0,	9,	0,	0,	0,	0,	0	</R>\n				<R>	Q,	0,	8,	2,	2,	0,	0,	5,	0,	0,	0,	0,	0	</R>\n				<R>	Y,	0,	9,	2,	1,	0,	0,	9,	0,	0,	0,	0,	0	</R>\n				<R>	M,	0,	9,	1,	1,	0,	1,	5,	0,	0,	0,	0,	0	</R>\n				<R>	H,	0,	7,	1,	3,	0,	0,	6,	0,	0,	0,	0,	0	</R>\n				<R>	C,	0,	5,	1,	1,	0,	1,	3,	0,	0,	0,	0,	0	</R>\n				<R>	W,	0,	10,	1,	2,	0,	0,	11,	0,	0,	0,	0,	0	</R>\n"

ptm_iso_free_prorata={"PTM{~} = M":"\t\t\t\t<R> ~,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{!} = NQR":"\t\t\t\t<R> !,  0,  -1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{@} = STYHD":"\t\t\t\t<R> @,  0,	1,	3,	0,	1,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{>to1} = STYHD":"\t\t\t\t<R> 1,  0,	0,	0,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{<to2} = ST":"\t\t\t\t<R> 2,  0,	-2,	-1,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{%} = K":"\t\t\t\t<R> %,  2,	2,	1,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{^} = KRED":"\t\t\t\t<R> ^,  1,	2,	0,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{&} = KR":"\t\t\t\t<R> &,  2,	4,	0,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{*} = K":"\t\t\t\t<R> *,  3,	6,	0,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{(} = C":"\t\t\t\t<R> (,  0,	-1,	1,	1,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{)} = Y":"\t\t\t\t<R> ),  0,	-1,	2,	1,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{/} = C":"\t\t\t\t<R> /,  2,	3,	1,	1,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{$} = D":"\t\t\t\t<R> $,  1,	2,	0,	0,	0,	1,  0,  0,  0,  0,  0,  0   </R>\n"}

ptm_iso_c_prorata={"PTM{~} = M":"\t\t\t\t<R> ~,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{!} = NQR":"\t\t\t\t<R> !,  0,  -1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{@} = STYHD":"\t\t\t\t<R> @,  0,	1,	3,	0,	1,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{>to1} = STYHD":"\t\t\t\t<R> 1,  0,	0,	0,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{<to2} = ST":"\t\t\t\t<R> 2,  0,	-2,	-1,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{%} = K":"\t\t\t\t<R> %,  0,	2,	1,	0,	0,	0,  2,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{^} = KRED":"\t\t\t\t<R> ^,  0,	2,	0,	0,	0,	0,  1,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{&} = KR":"\t\t\t\t<R> &,  0,	4,	0,	0,	0,	0,  2,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{*} = K":"\t\t\t\t<R> *,  0,	6,	0,	0,	0,	0,  3,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{(} = C":"\t\t\t\t<R> (,  0,	-1,	1,	1,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{)} = Y":"\t\t\t\t<R> ),  0,	-1,	2,	1,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{/} = C":"\t\t\t\t<R> /,  0,	3,	1,	1,	0,	0,  2,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{$} = D":"\t\t\t\t<R> $,  0,	2,	0,	0,	0,	1,  1,  0,  0,  0,  0,  0   </R>\n"}

ptm_iso_n_prorata={"PTM{~} = M":"\t\t\t\t<R> ~,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{!} = NQR":"\t\t\t\t<R> !,  0,  -1,  1, 0,  0,  0,  0,  0,  0,  -1,  0,  0   </R>\n",
                        "PTM{@} = STYHD":"\t\t\t\t<R> @,  0,	1,	3,	0,	1,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{>to1} = STYHD":"\t\t\t\t<R> 1,  0,	0,	0,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{<to2} = ST":"\t\t\t\t<R> 2,  0,	-2,	-1,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{%} = K":"\t\t\t\t<R> %,  2,	2,	1,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{^} = KRED":"\t\t\t\t<R> ^,  1,	2,	0,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{&} = KR":"\t\t\t\t<R> &,  2,	4,	0,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{*} = K":"\t\t\t\t<R> *,  3,	6,	0,	0,	0,	0,  0,  0,  0,  0,  0,  0   </R>\n",
                        "PTM{(} = C":"\t\t\t\t<R> (,  0,	-1,	1,	0,	0,	0,  0,  0,  0,  1,  0,  0   </R>\n",
                        "PTM{)} = Y":"\t\t\t\t<R> ),  0,	-1,	2,	0,	0,	0,  0,  0,  0,  1,  0,  0   </R>\n",
                        "PTM{/} = C":"\t\t\t\t<R> /,  2,	3,	1,	0,	0,	0,  0,  0,  0,  1,  0,  0   </R>\n",
                        "PTM{$} = D":"\t\t\t\t<R> $,  1,	2,	0,	0,	0,	1,  0,  0,  0,  0,  0,  0   </R>\n"}
