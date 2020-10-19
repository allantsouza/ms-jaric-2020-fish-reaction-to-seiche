# Input: raw hobo data
# Output: table of temperature by 0.1 m

hobo_data <- as_tibble(load_hobo_data())


temperatures_monotonic <- hobo_data %>% 
  group_nest(location, interval) %>%
  rename(ts = interval) %>%
  filter(map_lgl(data, ~ length(.x$depth) > 1 & all(!is.na(.x$depth & !is.na(.x$temperature))))) %>%
  mutate(data_hypo_extended = map(data, ~ as_tibble(approx(x = .x$depth, 
                                                         y  = .x$temperature,
                                                         xout = c(0, .x$depth, max(.x$depth) + 1:2),
                                                         rule = 2)))) %>%
  mutate(interpolated_profile = map(data_hypo_extended, ~ as_tibble(interpolate_temperature_profile(depth = .x$x, temperature = .x$y, depth_res = 0.1)))) %>%
  dplyr::select(-data, -data_hypo_extended) %>%
  unnest("interpolated_profile") %>%
  write_csv(file = here("data", "products", "temperature_data.csv"))
