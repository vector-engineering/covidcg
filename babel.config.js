module.exports = {
  plugins: [
    [
      '@babel/plugin-proposal-decorators',
      {
        legacy: true,
      },
    ],
    ['babel-plugin-styled-components'],
  ],
};
