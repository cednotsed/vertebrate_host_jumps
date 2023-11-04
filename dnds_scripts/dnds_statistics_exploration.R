pal <- distinctColorPalette(11)

test_df <- foreach(ds = seq(0, 0.5, 0.1), .combine = "bind_rows") %do% {
  tibble(dn = dn,
         ds = ds) %>%
    mutate(dnds = dn / ds)
}

plt1 <- test_df %>%
  ggplot(aes(x = dn, y = dnds, color = factor(ds))) +
  scale_color_manual(values = pal) +
  geom_line() +
  labs(y = "dN/dS")

test_df <- foreach(ds = seq(0, 0.5, 0.1), .combine = "bind_rows") %do% {
  tibble(dn = dn,
         ds = ds) %>%
    mutate(dnds = dn / (dn + ds))
}

plt2 <- test_df %>%
  ggplot(aes(x = dn, y = dnds, color = factor(ds))) +
  scale_color_manual(values = pal) +
  geom_line() +
  labs(y = "dN/(dN + dS)")

test_df <- foreach(ds = seq(0, 0.5, 0.1), .combine = "bind_rows") %do% {
  tibble(dn = dn,
         ds = ds) %>%
    mutate(dnds = dn / (ds + 1))
}

plt3 <- test_df %>%
  ggplot(aes(x = dn, y = dnds, color = factor(ds))) +
  scale_color_manual(values = pal) +
  geom_line() +
  labs(y = "dN/(dS + 1)")

test_df <- foreach(ds = seq(0, 0.5, 0.1), .combine = "bind_rows") %do% {
  tibble(dn = dn,
         ds = ds) %>%
    mutate(dnds = dn / (ds + 1e-6))
}

plt4 <- test_df %>%
  ggplot(aes(x = dn, y = dnds, color = factor(ds))) +
  scale_color_manual(values = pal) +
  geom_line() +
  labs(y = "dN/(dS + 1e-6)")


ggarrange(plt1, plt2, plt3, plt4, common.legend = T)

