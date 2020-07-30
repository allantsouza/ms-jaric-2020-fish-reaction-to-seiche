explain
                          WITH posth AS (SELECT th.*, lake, up_tag_sn, up_disttoshore, up_timestamp_utc FROM 
                          (SELECT up_tag_sn, up_id, up_disttoshore, up_timestamp_utc, up_depth, up_lake lake FROM at_macfish.umap_pos_powfilter WHERE
                              up_validpos AND fishvalid AND up_tag_sn IN ('T412150')) pos 
                          INNER JOIN 
                            (SELECT * FROM at_macfish.umap_pos_thermocline WHERE pth_step_order = 1 AND pth_degree_per_meter = 2 
                                                                                AND pth_thermpart IN ('center','start','end')
                              ) as th
                          ON pos.up_id = th.up_id )
                          
                          SELECT pos.*, dd_depth, dd_timestamp_utc, up_timestamp_utc FROM at_macfish.detsdepth dd 
                          INNER JOIN (SELECT posth.*, b.dd_id FROM posth INNER JOIN at_macfish.detsdepth_to_umap_pos b ON b.up_id = posth.up_id) pos ON
                          pos.dd_id = dd.dd_id;
