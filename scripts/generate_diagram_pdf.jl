root, dirs, files = (first∘collect∘walkdir)(".")

tex_list = filter( s->endswith(s,".tex"), files )
tex_head_list = map( s->s[1:end-4], tex_list )

pdf_list = filter( s->endswith(s,".pdf"), files )
pdf_head_list = map( s->s[1:end-4], pdf_list )

mission_head_list = setdiff( tex_head_list, pdf_head_list )

for one_mission in mission_head_list
  run( `lualatex $(one_mission)` )
end # for one_mission
