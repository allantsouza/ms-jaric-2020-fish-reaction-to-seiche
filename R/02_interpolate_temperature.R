# Input: raw hobo data
# Output: table of temperature by 0.1 m

hobo_data <- as_tibble(load_hobo_data())


temperatures_monotonic <- hobo_data %>% 
  group_nest(location, interval) %>%
  mutate(interpolated_profile = map(data, ~ as_tibble(interpolate_temperature_profile(depth = .x$depth, temperature = .x$temperature, depth_res = 0.1)))) %>%
  dplyr::select(-data) %>%
  unnest("interpolated_profile") %>%
  write_csv(file = here("data", "products", "temperature_data.csv"))
