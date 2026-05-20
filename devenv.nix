{
  pkgs,
  ...
}:

{
  packages = with pkgs; [
    entr
    fd
    ruff
    ty
  ];

  languages.python = {
    enable = true;
    package = pkgs.python3.withPackages (ps: [ ps.pytest ]);
    lsp.package = pkgs.ty;
    venv.enable = true;
    uv = {
      enable = true;
      sync = {
        enable = true;
        allGroups = true;
      };
    };
  };

  scripts.dev = {
    exec = ''
      fd -tf | entr -c pytest
    '';
    description = "run pytest on src file change";
    packages = with pkgs; [
      fd
      entr
    ];
  };

  enterShell = ''
    echo "========================================="
    echo "| Available commands:                   |"
    echo "|   dev - run pytest on src file change |"
    echo "========================================="
  '';

  # git-hooks.hooks.ruff.enable = true;
  # git-hooks.hooks.ruff-format.enable = true;
}
